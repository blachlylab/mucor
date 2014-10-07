#!/usr/bin/env python
# -*- coding: utf8
#
# mucor.py
# (c) James S Blachly, MD 2013
# 

# let print() be a function rather than statement
# ala python3
from __future__ import print_function

import os
import sys
import time
import argparse #transitioned from getopt
import csv
import itertools
import HTSeq
from collections import defaultdict
import gzip
import cPickle as pickle
import pdb   #pdb.set_trace()

import xml.etree.ElementTree as ET
import json

import numpy as np
import pandas as pd

# mucor modules
import mucorfilters as mf
from variant import Variant
from mucorfeature import MucorFeature
from databases import load_db, isAnnotatedSNP 

class Info:
    '''Program info: logo, version, and usage'''
    logo = """
 __    __     __  __     ______     ______     ______    
/\ "-./  \   /\ \/\ \   /\  ___\   /\  __ \   /\  == \   
\ \ \-./\ \  \ \ \_\ \  \ \ \____  \ \ \/\ \  \ \  __<   
 \ \_\ \ \_\  \ \_____\  \ \_____\  \ \_____\  \ \_\ \_\ 
  \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_/ /_/ 
                                  
"""
    version = "0.17"
    versionInfo = "mucor version {0}\nJames S Blachly, MD\nKarl W Kroll, BS".format(version)
    usage = """
Usage:
{0} [-h] | -g featurefile.gff -f feature_type [-u] -o <output_dir> <mutect001.txt mutect002.txt ... mutectNNN.txt>

Flags:
    -h  Show this help

    -g  GFF3/GTF file describing the features into which mutations will be binned

    -f  String to describe the type of feature to bin for this run.
        e.g. gene_id or transcript_id or chromosome_id

    -u,
    --union
        Join all items with same ID for feature_type (specified by -f)
        into a single, continuous bin. For example, if you want intronic
        variants counted with a gene, use this option. 
        ** TO DO **
        WARNING, this will likely lead to spurious results due to things
        like MIR4283-2 which exists twice on + and - strand of the same
        chromosome over 1 megabase apart. This creates one huge spurious
        bin.

    -o  Specify output directory

Output directory:
    Output files will be placed in the specified directory.
    If the directory already exists, an error will occur (won't overwrite)

    Output consists of CSV spreadsheets (.txt):
        1. Plain text report of # variants binned according to feature_type
        2. Summary of pre-filter and post-filter variants
        3. Detailed report of all variants by feature_type

Input files:
    <mutect001.txt mutect002.txt ... mutectNNN.txt>
    Final arguments should be a list of mutations in muTect output format

"""
    description = """

mucor: MUtation CORrelation

mucor reads in variant files from a variety of sources (VCF, muTect .out)
and counts the number of mutations falling into known features. These are
grouped together and output to see which features (genes) and which spec-
-ific locations within those genes have the highest frequency of mutation
within a group.
"""
    epilog = """

Output directory:
    Output files will be placed in the specified directory.
    If the directory already exists, an error will occur (won't overwrite)

    Output consists of CSV spreadsheets (.txt):
        1. Plain text report of # variants binned according to feature_type
        2. Summary of pre-filter and post-filter variants
        3. Detailed report of all variants by feature_type

Input files:
    <mutect001.txt mutect002.txt ... mutectNNN.txt>
    Final arguments should be a list of mutations in muTect output format
"""

def abortWithMessage(message):
    print("*** FATAL ERROR: " + message + " ***")
    exit(2)

def str_to_bool(s):
    if str(s) == 'False':
        return False
    else:
        return True

######## Karl Modified ##############
# new, separate function to construct the genomic array of sets
def constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures):
    gas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for feature in itertools.islice(gffFile, 0, None):
        # Nonstandard contigs (eg chr17_ctg5_hap1, chr19_gl000209_random, chrUn_...)
        # must be specifically excluded, otherwise you will end up with exception
        # ValueError: start is larger than end due to duplicated gene symbols
        if "_hap" in feature.iv.chrom:
            continue
        elif "_random" in feature.iv.chrom:
            continue
        elif "chrUn_" in feature.iv.chrom:
            continue

        # transform feature to instance of Class MucorFeature
        feat = MucorFeature(feature.attr[featureType], feature.type, feature.iv)
        
        # WARNING
        # the following REQUIRES a coordinate-sorted GFF/GTF file
        # extra checks incurring slowdown penalty are req'd if GFF/GTF not sorted
        if feat.name in knownFeatures:
            # In case there is an error in the GFF and/or the featureType (-f) is not unique,
            # issue a warning
            # for example, genes.gtf supplied with the Illumina igenomes package for the tuxedo tools suite
            # includes duplicate entries for many genes -- e.g. DDX11L1 on chr15 shoudl be DDX11L9
            # try to cope with this by relabeling subsequent genes as GENESYM.chrNN
            if feat.iv.chrom != knownFeatures[feat.name].iv.chrom:
                duplicateFeatures.add(feat.name)
                feat.name = feat.name + '.' + feat.iv.chrom
            else:
                # do not obliterate the start coordinate when adding SUCCESSIVE bits of a feature (e.g. exons)
                # TO DO: Here is where the --nounion option would work
                #feat.iv.start = knownFeatures[feat.name].iv.start
                pass # no-union - this does overwrite previous coordinates in knownFeatures,
                     # but should not matter as the actual coordinates are obtaind from 'gas'.

        # first, add to the knownFeatures, a dictionary of MucorFeatures, which contain the variants set
        knownFeatures[feat.name] = feat
        # then, add to the GenomicArrayOfSets, which we use to find gene symbol from variant coords
        try:
            gas[ feat.iv ] += feat.name
        except ValueError:
            print(feat.name)
            print(feat.iv)
            raise
    return gas, knownFeatures, duplicateFeatures

def parseJSON(json_config):
    global filename2samples
    filename2samples = {}
    JD = json.load(open(json_config,'r'))
    featureType = JD['feature']
    outputDir = JD['run_name']
    union = JD['union']
    fast = JD['fast']
    gff = JD['gtf']
    database = JD['database']
    filters = JD['filters']
    global database_switch
    database_switch = str_to_bool(database[0])
    global SnpEff_switch
    SnpEff_switch = False
    '''
    for i in JD['filters']:
        filters[i] = True      # Imagine filters as "ON/OFF", binary switches
    '''
    input_files = []
    for i in JD['samples']:
        for j in i['files']:
            filename = str(j['path']).split('/')[-1]
            filename2samples[filename] = i['id']
            '''
            if not MuTect_switch and str(j['source']) == str('Mutect') and str(j['type']) == str('mutect'):
                MuTect_switch = bool(True)
                pass
            '''
            if not SnpEff_switch and str(j['type']) == str('vcf') and bool(j['snpeff']) == bool(True):
                SnpEff_switch = bool(True)
            input_files.append(j['path'])

    return featureType, outputDir, union, fast, gff, database, filters, input_files

def parseGffFile(gffFileName, featureType, fast):
    '''Parse the GFF/GTF file. Return tuple (knownFeatures, GenomicArrayOfSets)
    Haplotype contigs are explicitly excluded because of a coordinate crash (begin > end)'''
    
    # TO DO: command line flag should indicate that variants in INTRONS are counted
    # This is called --union, see below
    
    startTime = time.clock()
    print("\n=== Reading GFF/GTF file {0} ===".format(gffFileName))
    
    scripts = "/".join(os.path.realpath(__file__).split('/')[:-1])
    annotFileName = ".".join(gffFileName.split('/')[-1].split('.')[:-1])
    archiveFilePath = str("/") + str(scripts).strip('/') + str("/") + str(annotFileName) + str('_') + str(featureType) + str('.p')
    print(gffFileName)
    gffFile = HTSeq.GFF_Reader(gffFileName)
    #ga = HTSeq.GenomicArray("auto", typecode="i")  # typecode i is integer
    
    knownFeatures = {}                              # empty dict

    duplicateFeatures = set()

    # gas - GenomicArrayOfSets.
    # typecode always 'O' (object) for GenomicArrayOfSets
    # UNstranded -- VCF and muTect output always report on + strand,
    # but the GenomicArray must be unstranded because the GFF /is/ strand-specific,
    # and if I manually coded all GenomicIntervals read from the VCF or muTect file as '+',
    # then no genes on the - strand would have variants binned to them

    if fast:
        if os.path.exists( archiveFilePath ):
            print("Opening annotation archive: " + str(archiveFilePath))
            pAnnot = open(archiveFilePath, 'rb')
            gas = pickle.load(pAnnot)
            knownFeatures = pickle.load(pAnnot)
            duplicateFeatures = pickle.load(pAnnot)
            pAnnot.close()
        else:
            print("Cannot locate annotation archive for " + str(gffFileName.split('/')[-1]) + str(" w/ ") + str(featureType))
            print("   Reading in annotation and saving archive for faster future runs") 
            gas, knownFeatures, duplicateFeatures = constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures)
            archiveOut = open(archiveFilePath, 'wb')
            pickle.dump(gas, archiveOut, -1) ### pickle feature only works with full annotation files
            pickle.dump(knownFeatures, archiveOut, -1)
            pickle.dump(duplicateFeatures, archiveOut, -1)
            archiveOut.close()
    if not fast:
    # ignore pickles function
        gas, knownFeatures, duplicateFeatures = constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures)

    if duplicateFeatures:
        print("*** WARNING: {0} {1}s found on more than one contig".format(len(duplicateFeatures), featureType))

    totalTime = time.clock() - startTime
    print("{0} sec\t{1} found:\t{2}".format(int(totalTime), featureType, len(knownFeatures)))

    return knownFeatures, gas

def parse_MiSeq(row, fieldId, header):
    VF = row[fieldId[header[-1]]].split(':')[-2]
    DP = row[fieldId[header[-1]]].split(':')[2]
    position = int(row[fieldId['POS']])
    return VF, DP, position

def parse_IonTorrent(row, fieldId, header):
    for i in row[fieldId['INFO']].split(';'):
        if i.startswith("AO="):
            tempval = i.split('=')[1]
        if i.startswith("RO="):
            RO = i.split('=')[1]
        if i.startswith("DP="):
            DP = i.split("=")[1]
    if str(',') in str(tempval):
        tempval2 = [int(numeric_string) for numeric_string in tempval.split(',')]
        try:
            AO = sum(tempval2)
        except:
            abortWithMessage("AO should be an int, or a list of ints: AO = {0}/".format(tempval2))
    else:
        AO = tempval
    VF = float(float(AO)/float(float(RO) + float(AO)))
    position = int(row[fieldId['POS']])
    for i in str(row[fieldId['ALT']]).split(','):
        if len(str(row[fieldId['REF']])) > len(i):
            #this is a deletion in Ion Torrent data
            position = int(row[fieldId['POS']])
            break
    return VF, DP, position

def parse_MuTectOUT(row, fieldId): # ,MuTect_Annotations
    VF = row[fieldId['tumor_f']]
    DP = int(int(str(row[fieldId['t_ref_count']]).strip()) + int(str(row[fieldId['t_alt_count']]).strip()))
    position = int(row[fieldId['position']])
    '''
    MuTect_Annotations[tuple(( str(row[0]), position) )] = row[fieldId['dbsnp_site']]
    '''

    return VF, DP, position # , MuTect_Annotations 
    
def parse_MuTectVCF(row, fieldId, header, fn): # , MuTect_Annotations)
    j = 0
    for i in header:
        if str('-') in str(i): # This line should detect if the sample id is in the line. Should be rewritten for samples without a "-" in their name
            tmpsampID = i
    for i in row[fieldId['FORMAT']].split(':'):
        if i == "FA":
            VF = row[fieldId[tmpsampID]].split(':')[j]
        elif i == "DP":
            DP = row[fieldId[tmpsampID]].split(':')[j]
        j+=1
    position = int(row[fieldId['POS']])
    '''
    global MuTect_switch
    MuTect_switch = True
    
    if MuTect_switch == True and os.path.exists(fn.replace('_snpEff.vcf', '.out')) and str(fn.replace('_snpEff.vcf', '.out')) != str(fn):
        MuTect_output = fn.replace('_snpEff.vcf', '.out')
        for line in open(MuTect_output):
            if str(str(row[0]) + "\t") in str(line) and str(str(position) + "\t") in str(line):
                MuTect_Annotations[tuple((str(row[0]), position))] = line.split('\t')[8]
                break
            else:
                continue
    '''

    return VF, DP, position # , MuTect_Annotations

def parse_SomaticIndelDetector(row, fieldId, header):
    j = 0
    for i in header:
        if str('-') in str(i):
            tmpsampID = i
    for i in row[fieldId['FORMAT']].split(':'):
        if i == "AD":
            ALT_count = row[fieldId[tmpsampID]].split(':')[j].split(',')[1]
        elif i == "DP":
            DP = row[fieldId[tmpsampID]].split(':')[j]
            VF = float( float(ALT_count)/float(DP) )
        j+=1
    position = int(row[fieldId['POS']])
    return VF, DP, position

def parse_SamTools(row, fieldId, header):
    position = int(row[fieldId['POS']])
    for i in row[fieldId['INFO']].split(';'):
        if i.startswith("DP4="):
            j = i.split('=')[1].split(',')
            ref = int(int(j[0]) + int(j[1]))
            alt = int(int(j[2]) + int(j[3]))
            DP = int(int(ref) + int(alt))
            VF = float( float(alt)/float(DP) )
    return VF, DP, position

def parse_VarScan(row, fieldId, header):
    j = 0
    position = int(row[fieldId['POS']])
    for i in row[fieldId['FORMAT']].split(':'):
        if str(i) == "DP":
            DP = int(row[fieldId[header[-1]]].split(':')[j])
        if str(i) == "FREQ":
            VF = float(float(str(row[fieldId[header[-1]]].split(':')[j]).strip('%'))/float(100))
        j += 1
    return VF, DP, position

def filterRow(row, fieldId, filters):
    try:
        for rowFilter in str(row[fieldId['FILTER']]).split(';'):    ## VCF file format
            if rowFilter not in filters:
                return True
                break
    except KeyError:
        for rowFilter in str(row[fieldId['judgement']]).split(';'): ## MuTect '.out' file format
            if rowFilter not in filters:
                return True
                break
    return False

def parseVariantFiles(variantFiles, knownFeatures, gas, snps, filters):
    # parse the variant files (VCF, muTect format)
    startTime = time.clock()

    # All variants stored in long (record) format
    # in a pandas dataframe
    varD = {'chr':[],'pos':[],'ref':[],'alt':[],'vf':[],'dp':[],'feature':[],'effect':[],'fc':[],'datab':[],'sample':[],'source':[]}

    print("\n=== Reading Variant Files ===")
    for fn in variantFiles:

        varFile = open(fn, 'rb')    # TO DO: error handling
        varReader = csv.reader(varFile, delimiter='\t')

        try:
            row = varReader.next()
        except StopIteration:
            # a file was empty (i.e. no first row to read)
            print("Empty file {}".format(fn))
            continue    # next fn in variantFiles
        #if len(row) != 1: raise ValueError('Invalid muTector header')
        #if "## muTector" not in row[0]: raise ValueError('Invalid muTector header')
        global MiSeq
        MiSeq = False
        global IonTorrent
        IonTorrent = False
        global Mutect
        Mutect = False
        global Mutector
        Mutector = False
        global SomaticIndelDetector
        SomaticIndelDetector = False
        global Samtools
        Samtools = False
        global VarScan 
        VarScan = False
        while str(row).split("'")[1][0:2] == '##':
            if str('Torrent Unified Variant Caller') in str(row): 
                IonTorrent = True
            elif str('MiSeq') in str(row):
                MiSeq = True
            elif str('SomaticIndelDetector') in str(row):
                SomaticIndelDetector = True
            elif str('MuTect') in str(row):
                Mutect = True
            elif str('muTector') in str(row):
                Mutector = True
            elif str('samtools') in str(row):
                Samtools = True
            elif str('source=VarScan') in str(row):
                VarScan = True
            row = varReader.next()
        
        header = row
        if len(header) == 0: raise ValueError('Invalid header')
        fieldId = dict(zip(header, range(0, len(header))))

        # read coverage depth minimum cutoff; currently unusued
        read_depth_min = 0
        # after reading the two header rows, read data
        for row in itertools.islice(varReader, None):
            if filterRow(row, fieldId, filters):    # filter rows as they come in, to prevent them from entering the dataframe
                continue                            # this allows us to print the dataframe directly and have consistent output with variant_details.txt, etc.

            # make a variant object for row
            # TO DO: change row index#s to column names or transition to row object

            EFF = ""
            FC = ""
            muts = []
            loca = []
            try:
                for eff in row[fieldId['INFO']].split(';'):
                    if eff.startswith('EFF='):
                        for j in eff.split(','):
                            muts.append(str(j.split('|')[3]))
                            loca.append(str(j.split('(')[0]).replace('EFF=',''))
                for guy in set(muts):
                    if str(guy) != "":
                        EFF += str(guy) + ";"
                for guy in set(loca):
                    if str(guy) != "":
                        FC += str(guy) + ";"
            except KeyError:
                pass
            if MiSeq:
                VF, DP, position = parse_MiSeq(row, fieldId, header)
            elif IonTorrent:
                VF, DP, position = parse_IonTorrent(row, fieldId, header)
            elif Mutect:
                VF, DP, position = parse_MuTectVCF(row, fieldId, header, fn) # , MuTect_Annotations)
            elif SomaticIndelDetector:
                VF, DP, position = parse_SomaticIndelDetector(row, fieldId, header)
            elif Mutector:
                VF, DP, position = parse_MuTectOUT(row, fieldId) # , MuTect_Annotations)
            elif Samtools:
                VF, DP, position = parse_SamTools(row, fieldId, header)
            elif VarScan:
                VF, DP, position = parse_VarScan(row, fieldId, header)
            else:
                print("This isn't MiSeq, IonTorrent, SomaticIndelDetector, Samtools, VarScan, or Mutect data?")
                sys.exit(1)
            var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(row[0], int(position)), ref=row[3], alt=row[4], frac=VF, dp=DP, eff=EFF.strip(';'), fc=FC.strip(';'))
            ###########################################

            # find bin for variant location
            resultSet = gas[ var.pos ]      # returns a set of zero to n IDs (e.g. gene symbols)
            if resultSet:                   # which I'll use as a key on the knownFeatures dict
                #print var.pos              # and each feature with matching ID gets allocated the variant
                #print(gas[ var.pos ])      # 
                for featureName in resultSet:
                    knownFeatures[featureName].variants.add(var)
            
            # Descriptive variable names
            chr = row[0]
            pos = int(position)
            ref = row[3]
            alt = row[4]
            vf = VF     # to do, change (karl capitalizes)
            dp = DP     # to do, change (karl capitalizes)
            feature = ', '.join( gas[ var.pos ] )   # join with comma to handle overlapping features
            effect = EFF
            fc = FC
            #sample = fn.split('/')[-1]      # to do, will need to come from JSON config
            sample = filename2samples[str(fn.split('/')[-1])]
            source = fn.split('/')[-1]
            if snps.has_key((chr, pos)):
                datab = str( x for x in set(snps[(chr,pos)][0]) )
            else:
                datab = str('?')
            '''
            if Mutect or Mutector:
                datab = MuTect_Annotations[(chr, pos)]
            '''
            # build dict to insert
            #columns=('chr','pos','ref','alt','vf','dp','gene','effect','sample','source')
            vardata = dict(zip( ['chr','pos','ref','alt','vf','dp','feature','effect','fc','datab','sample','source'], \
                                [ chr , pos , ref , alt , vf , dp , feature , effect , fc , datab , sample , source ])) 
            # add to variants data frame
            for key in vardata.keys():
                varD[key].append(vardata[key])

        totalTime = time.clock() - startTime
        print("{0:02d}:{1:02d}\t{2}".format(int(totalTime/60), int(totalTime % 60), fn))
    
    # Removing columns from the following 'columns' list will mask them from output in allvars.txt
    varDF = pd.DataFrame(varD, columns=['chr','pos','ref','alt','vf','dp','feature','effect','fc','datab','sample','source'])
    # Clean up variant dataframe a little
    # position should be integer, not float
    varDF.pos = varDF.pos.astype(int)
    return varDF, knownFeatures, gas, snps

def printOutput(argv, outputDirName, knownFeatures, gas, snps): ######## Karl Modified ##############
    '''Output statistics and variant details to the specified output directory.'''

    startTime = time.clock()
    print("\n=== Writing output files to {0}/ ===".format(outputDirName))

    try:
        # of = outputFile
        ofRunInfo = open(outputDirName + "/run_info.txt", 'w+')
        ofCounts = open(outputDirName + "/counts.txt", 'w+')
        ofVariantDetails = open(outputDirName + "/variant_details.txt", 'w+')
        ofVariantBeds = open(outputDirName + "/variant_locations.bed", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))

    # =========================
    # run_info.txt
    #
    ofRunInfo.write(Info.versionInfo + '\n')
    ofRunInfo.write("{0}\n".format(time.ctime() ) )
    ofRunInfo.write("Command line: {0}\n".format(str(argv) ) )
    ofRunInfo.write("No. samples: \n")
    ofRunInfo.write("Filters:\n")
    #
    ofRunInfo.write("Variants Pre-filter: \n")
    ofRunInfo.write("        Post-filter: \n")
    ofRunInfo.close()

    # ============================================================
    # counts.txt
    # make master list, then sort it by number of variants per bin
    #
    ofCounts.write('FeatureName\tHits\tWeightedHits\tAverageWeight\tUniqueHits\tNumSamples\n')
    masterList = list(knownFeatures.values())
    sortedList = sorted(masterList, key=lambda k: k.numVariants(), reverse=True)
    nrow = 0

    for feature in sortedList:
        if knownFeatures[feature.name].variants:
            ofCounts.write(feature.name + '\t')

            ofCounts.write(str(len(knownFeatures[feature.name].variants)) + '\t')
            
            ofCounts.write(str(knownFeatures[feature.name].weightedVariants()) + '\t')
            
            ft = knownFeatures[feature.name]
            avgWt = float(ft.weightedVariants() / float(ft.numVariants()) )
            ofCounts.write(str(avgWt) + '\t')

            ofCounts.write(str(knownFeatures[feature.name].numUniqueVariants()) + '\t')
            
            ofCounts.write(str(knownFeatures[feature.name].numUniqueSamples()) + '\n')

            nrow += 1
        '''
        else:
            print(feature)
        '''
    
    print("\t{0}: {1} rows".format(ofCounts.name, nrow))
    ofCounts.close()
    totalTime = time.clock() - startTime
    print("\tTime to write: {0:02d}:{1:02d}".format(int(totalTime/60), int(totalTime % 60)))

    # =========================================================
    # variant_details.txt
    #
    ofVariantDetails.write('Feature\tContig\tPos\tRef\tAlt\tVF\tDP\tSource\t') ######## Karl Modified ##############
    if SnpEff_switch:
        ofVariantDetails.write('Effect\tFC\t')
    if database_switch:
        ofVariantDetails.write('Annotation\t')
    ofVariantDetails.write('Count\n')

    for feature in sortedList:
        if knownFeatures[feature.name].variants:
            for var in knownFeatures[feature.name].uniqueVariants():
                
                #if not isAnnotatedSNP(snps, tuple((var.pos.chrom,var.pos.pos))):  ##### KARL ADDED ######
                ofVariantDetails.write(feature.name + '\t')
                ofVariantDetails.write(var.pos.chrom + '\t')
                ofVariantBeds.write(var.pos.chrom + '\t')
                ofVariantDetails.write(str(var.pos.pos) + '\t')
                ofVariantBeds.write(str(var.pos.pos - 1) + '\t')
                ofVariantBeds.write(str(var.pos.pos) + '\t')
                #ofVariantBeds.write(str(snps[tuple((var.pos.chrom,var.pos.pos))]) + '\t')
                ofVariantDetails.write(var.ref + '\t')
                ofVariantDetails.write(var.alt + '\t')
                ofVariantDetails.write(var.frac + '\t')
                ofVariantDetails.write(var.dp + '\t')
                ofVariantDetails.write(var.source + '\t')
                if SnpEff_switch:
                    if str(var.eff) != str(''):
                        ofVariantDetails.write(str([ x for x in set(str(var.eff).split(', '))]).replace("''","").replace(", ","").strip(']').strip('[').strip("'") + '\t')
                    else:
                        ofVariantDetails.write(str('?') + '\t')
                    if str(var.fc) != str(''):
                        ofVariantDetails.write(str([ x for x in set(str(var.fc).split(', '))]).replace("''","").replace(", ","").strip(']').strip('[').strip("'") + '\t')
                    else:
                        ofVariantDetails.write(str('?') + '\t')
                if database_switch:
                    if isAnnotatedSNP(snps, tuple((str(var.pos.chrom),str(var.pos.pos)))):
                        ofVariantDetails.write(snps[((str(var.pos.chrom),str(var.pos.pos)))][0] + '\t') # snps[(var.pos.chrom, var.pos.pos)].rs
                    else:
                        ofVariantDetails.write(str('?') + '\t')
                ofVariantDetails.write(str(len(var.source.split(','))) + '\n')
                ofVariantBeds.write('\n')
        ########### Karl added bed file output here #############
    ofVariantDetails.close()
    ofVariantBeds.close()
    
    return 0


def main():
    print(Info.logo)
    print(Info.versionInfo)

    print("\n=== Run Info ===")
    print("\t{0}".format(time.ctime() ) )
    print()

    featureType, outputDir, union, fast, gff, database, filters, input_files = parseJSON(sys.argv[1])

    if not os.path.exists(gff):
        abortWithMessage("Could not find GFF file {0}".format(gff))
    if os.path.exists(outputDir) and os.listdir(outputDir):
        abortWithMessage("The directory {0} already exists and contains output. Will not overwrite.".format(outputDir))
    elif not os.path.exists(outputDir):
        try:
            os.makedirs(outputDir)
        except:
            abortWithMessage("Error when creating output directory {0}".format(outputDir))

    # check that all specified variant files exist
    for fn in input_files:
        if not os.path.exists(fn):
            abortWithMessage("Could not find variant file: {0}".format(fn))


    knownFeatures, gas = parseGffFile(str(gff), str(featureType), bool(fast))

    snps =  load_db(database) ######## Karl Modified ##############
    '''
    outsnps = open('outsnps.p','rb')
    pickle.dump(snps,outsnps)
    outsnps.close()
    abortWithMessage("done making snps dict")
    '''
    varDF, knownFeatures, gas, snps = parseVariantFiles(list(input_files), knownFeatures, gas, snps, filters)
    printOutput(list(input_files), str(outputDir), knownFeatures, gas, snps) ######## Karl Modified ##############
    
    ## ## ##
    # print record format (long format) all variants data frame
    # TO DO : break out into function to slicing and dicing
    #if (varDF['chr'][0].lower().startswith(('chr', 'Chr'))):
    #    varDF['chr_num'] = varDF['chr'].apply(lambda x:x[3:]).astype(int)   # strip off the 'chr', unstringify
    #    varDF.sort(columns=['chr_num', 'pos'], inplace=True)
    #   varDF.drop(['chr_num'], axis=1, inplace=True)
    #else: 
    
    # sorting by chr does not sort properly: lexical order is chr1,chr10,...,chr2,...
    # so we will sort by feature (gene etc.) then position
    varDF.sort(columns=['feature','pos'], inplace=True)
    varDF.replace('', np.nan, inplace=True)
    varDF.to_csv(outputDir + '/allvars.txt', sep='\t', na_rep='?', index=False)
    
    # pretty print newline before exit
    print()


if __name__ == "__main__":
    if sys.hexversion < 0x02070000:
        raise RuntimeWarning("mucor should be run on python 2.7.0 or greater.")
    main()
