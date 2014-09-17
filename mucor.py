#!/usr/bin/env python
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


######## Karl Added ##############
# makes a dictionary out of dbSNP,
# with tuple of chrom,position as the key; and the rs number as values.
#
def load_dbsnp():
    startTime = time.clock()
    snps = defaultdict(str)
    dbsnp_p = '/nfs/17/osu7366/projects/new_AK/dbSNPandMiSeq.P'
    #dbsnp_file = '/nfs/17/osu7366/reference/snp138Common.txt.gz'
    #dbsnp = gzip.open(dbsnp_file,'rb')
    print("\n=== Reading dbSNP pickle file {0} ===".format(dbsnp_p))
    '''
    for line in dbsnp:
        col = line.split('\t')
        if str(col[11]) == "deletion": # deletions in our VCF file start 1 base upstream (-1) from dbSNP, but have the correct rs number
            snps[tuple((str(col[1]), int(col[3]) - 1))] = str(col[4])
        else:
            snps[tuple((str(col[1]), int(col[3])))] = str(col[4])
    '''
    snps = pickle.load(open(dbsnp_p,'rb'))
    totalTime = time.clock() - startTime
    print("{0} sec\t{1} SNPs".format(int(totalTime), len(snps.values())))
    return snps

##################################

######## Karl Added ##############
# true or false to check if a location
# (tuple of chrom,position) is in the dbSNP dictionary. 
# must use defaultdict above to avoid key errors here
def in_dbsnp(snps, loc):
    status = False
    annotation = snps[loc]
    if str(annotation).startswith('rs'):
        status = True
    return status

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

def parseGffFile(gffFileName, featureType, fast):
    '''Parse the GFF/GTF file. Return tuple (knownFeatures, GenomicArrayOfSets)
    Haplotype contigs are explicitly excluded because of a coordinate crash (begin > end)'''
    
    # TO DO: command line flag should indicate that variants in INTRONS are counted
    # This is called --union, see below
    
    startTime = time.clock()
    print("\n=== Reading GFF/GTF file {0} ===".format(gffFileName))
    
    scripts = "/".join(os.path.realpath(__file__).split('/')[:-1])
    annotFileName = ".".join(gffFileName.split('/')[-1].split('.')[:-1])
    archiveFilePath = str("/") + str(scripts).strip('/') + str("/") + str(annotFileName) + str('.p')
    
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
            print("Cannot locate annotation archive for " + str(gffFileName.split('/')[-1]))
            print("   Reading in annotation and saving archive for faster future runs") 
            gas, knownFeatures, duplicateFeatures = constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures)
            archiveOut = open(archiveFilePath, 'wb')
            pickle.dump(gas, archiveOut) ### pickle feature only works with full annotation files
            pickle.dump(knownFeatures, archiveOut)
            pickle.dump(duplicateFeatures, archiveOut)
            archiveOut.close()
    if not fast:
	# ignore pickles function
        gas, knownFeatures, duplicateFeatures = constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures)

    if duplicateFeatures:
        print("*** WARNING: {0} {1}'s found on more than one contig".format(len(duplicateFeatures), featureType))

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
            print("what's up with this? " + str(tempval2) )
            sys.exit(1)
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

def parse_MuTect(row, fieldId, header, fn, MuTect_Annotations):
    j = 0
    for i in header:
        if str('-') in str(i):
            tmpsampID = i
    for i in row[fieldId['FORMAT']].split(':'):
        if i == "FA":
            VF = row[fieldId[tmpsampID]].split(':')[j]
        elif i == "DP":
            DP = row[fieldId[tmpsampID]].split(':')[j]
        j+=1
    position = int(row[fieldId['POS']])

    global MuTect_switch
    MuTect_switch = True
    if MuTect_switch == True:
        MuTect_output = fn.replace('_snpEff.vcf', '.out')
        for line in open(MuTect_output):
            if str(str(row[0]) + "\t") in str(line) and str(str(position) + "\t") in str(line):
                MuTect_Annotations[tuple((str(row[0]), position))] = line.split('\t')[8]
                break
            else:
                continue
    return VF, DP, position, MuTect_Annotations

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


def parseVariantFiles(variantFiles, knownFeatures, gas, snps):
    # parse the variant files (muTect format)
    # TO DO: also interpret from VCF

    startTime = time.clock()
    global MuTect_Annotations
    MuTect_Annotations = defaultdict(str)

    # All variants stored in long (record) format
    # in a pandas dataframe
    varD = {'chr':[],'pos':[],'ref':[],'alt':[],'vf':[],'dp':[],'feature':[],'effect':[],'datab':[],'sample':[],'source':[]}

    print("\n=== Reading Variant Files ===")
    for fn in variantFiles:
        #print("\t{0}\t".format(fn), end='')    # moved to end to print after clock

        varFile = open(fn, 'rb')    # TO DO: error handling
        varReader = csv.reader(varFile, delimiter='\t')

        # '## muTector v1.0.47986'
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
        global SomaticIndelDetector
        SomaticIndelDetector = False
        while str(row).split("'")[1][0:2] == '##':
            if str('Torrent Unified Variant Caller') in str(row): 
                IonTorrent = True
            elif str('MiSeq') in str(row):
                MiSeq = True
            elif str('SomaticIndelDetector') in str(row):
                SomaticIndelDetector = True
            elif str('MuTect') in str(row):
                Mutect = True
            row = varReader.next()
        
        header = row
        if len(header) == 0: raise ValueError('Invalid header')
        fieldId = dict(zip(header, range(0, len(header))))

        # read coverage depth minimum cutoff; currently unusued
        read_depth_min = 0
        # after reading the two header rows, read data
        for row in itertools.islice(varReader, None):

            ############## Karl Added / Modified ################
            ######### Added variant frequency and depth to the file reading ########

            #############
            ## FILTERS ##   
            #############
            #if row[fieldId['FILTER']] != 'PASS': continue
            if row[fieldId['FILTER']] == 'REJECT': continue
            '''
            if MiSeq and str(row[fieldId['ID']])[0:2] == 'rs': 
                chrom = str(row[fieldId['#CHROM']])
                position = str(row[fieldId['POS']])
                print("found annotated mutation " + str(row[fieldId['ID']]) + " not in snp dictionary\n\tadding it now")
                snps[tuple((str(chrom),int(position)))] = str(row[fieldId['ID']])
                continue
            '''
            #if str(row[fieldId['INFO']]).split(';')[3].split('=')[1] >= int(read_depth_min): continue
            
            # make a variant object for row
            # TO DO: change row index#s to column names or transition to row object

            EFF = ""
            muts = []
            for eff in row[fieldId['INFO']].split(';'):
                if eff.startswith('EFF='):
                    for j in eff.split(','):
                        muts.append(str(j.split('|')[3]))
            for guy in set(muts):
                if str(guy) != "":
                    EFF += str(guy) + ";"

            if MiSeq:
                VF, DP, position = parse_MiSeq(row, fieldId, header)
            elif IonTorrent:
                VF, DP, position = parse_IonTorrent(row, fieldId, header)
            elif Mutect:
                VF, DP, position, MuTect_Annotations = parse_MuTect(row, fieldId, header, fn, MuTect_Annotations)
            elif SomaticIndelDetector:
                VF, DP, position = parse_SomaticIndelDetector(row, fieldId, header)

            else:
                print("This isn't MiSeq, IonTorrent, SomaticIndelDetector, or Mutect data?")
                sys.exit(1)
            var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(row[0], int(position)), ref=row[3], alt=row[4], frac=VF, dp=DP, eff=EFF.strip(';'))
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
            sample = fn.split('/')[-1]      # to do, will need to come from JSON config
            source = fn.split('/')[-1]
            datab = str('?')
            if Mutect:
                datab = MuTect_Annotations[(chr, pos)]
            # build dict to insert
            #columns=('chr','pos','ref','alt','vf','dp','gene','effect','sample','source')
            vardata = dict(zip( ['chr','pos','ref','alt','vf','dp','feature','effect','datab','sample','source'], \
                                [ chr , pos , ref , alt , vf , dp , feature , effect , datab , sample , source ])) 
            # add to variants data frame
            for key in vardata.keys():
                varD[key].append(vardata[key])

        totalTime = time.clock() - startTime
        print("{0:02d}:{1:02d}\t{2}".format(int(totalTime/60), int(totalTime % 60), fn))
    
    # Clean up variant dataframe a little
    # position should be integer, not float
    varDF = pd.DataFrame(varD, columns=['chr','pos','ref','alt','vf','dp','feature','effect','datab','sample','source'])
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
    ofVariantDetails.write('Feature\tContig\tPos\tRef\tAlt\tVF\tDP\tEffect\tSource\t') ######## Karl Modified ##############
    if MuTect_switch:
        ofVariantDetails.write('Annotation\tCount\n')
    elif not MuTect_switch:
        ofVariantDetails.write('Count\n')
    for feature in sortedList:
        if knownFeatures[feature.name].variants:
            for var in knownFeatures[feature.name].uniqueVariants():
                if not in_dbsnp(snps, tuple((var.pos.chrom,var.pos.pos))):  ##### KARL ADDED ######
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
                    ofVariantDetails.write(str([ x for x in set(str(var.eff).split(', '))]).strip(']').strip('[').strip("'") + '\t')
                    ofVariantDetails.write(var.source + '\t')
                    if MuTect_switch:
                        ofVariantDetails.write(MuTect_Annotations[(var.pos.chrom, var.pos.pos)] + '\t')
                    ofVariantDetails.write(str(len(var.source.split(','))) + '\n')
                    ofVariantBeds.write('\n')
        ########### Karl added bed file output here #############
    ofVariantDetails.close()
    ofVariantBeds.close()
    
    return 0


def main(argv):
    print(Info.logo)
    print(Info.versionInfo)

    print("\n=== Run Info ===")
    print("\t{0}".format(time.ctime() ) )
    print()

    parser = argparse.ArgumentParser(description=Info.description, epilog=Info.epilog, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-g", "--gff", required=True, help="Annotation GFF/GTF for feature binning")
    parser.add_argument("-f", "--featuretype", required=True, help="Feature type into which to bin [gene]")
    parser.add_argument("-u", "--union", action="store_true", help="""
        Join all items with same ID for feature_type (specified by -f)
        into a single, continuous bin. For example, if you want intronic
        variants counted with a gene, use this option. 
        ** TO DO **
        WARNING, this will likely lead to spurious results due to things
        like MIR4283-2 which exists twice on + and - strand of the same
        chromosome over 1 megabase apart. This creates one huge spurious
        bin.
        """)
    
    parser.add_argument("-o", "--output", required=True, help="Output directory name")
    parser.add_argument("variantfiles", nargs='+', help="List of variant files to muCorrelate")
    parser.add_argument("-n", "--no_archive", action="store_false", default=True, help="prevent quick load of annotation files")
    
    args = parser.parse_args()

    if not os.path.exists(args.gff):
        abortWithMessage("Could not find GFF file {0}".format(args.gff))
    if os.path.exists(args.output):
        abortWithMessage("The directory {0} already exists. Will not overwrite.".format(args.output))
    else:
        try:
            os.makedirs(args.output)
        except:
            abortWithMessage("Error when creating output directory {0}".format(outputDirName))

    # check that all specified variant files exist
    for fn in args.variantfiles:
        if not os.path.exists(fn):
            abortWithMessage("Could not find variant file: {0}".format(fn))

    fast = args.no_archive

    knownFeatures, gas = parseGffFile(args.gff, args.featuretype, fast)

    snps =  defaultdict(tuple) # load_dbsnp() ######## Karl Added ##############

    varDF, knownFeatures, gas, snps = parseVariantFiles(args.variantfiles, knownFeatures, gas, snps)
    
    printOutput(argv, args.output, knownFeatures, gas, snps) ######## Karl Modified ##############
    
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
    varDF.to_csv(args.output + '/allvars.txt', sep='\t', na_rep='?', index=False)
    
    # pretty print newline before exit
    print()


if __name__ == "__main__":
    if sys.hexversion < 0x02070000:
        raise RuntimeWarning("mucor should be run on python 2.7.0 or greater.")
    main(sys.argv[1:])
