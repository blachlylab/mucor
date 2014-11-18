#!/usr/bin/env python
# -*- coding: utf8
#
# mucor.py
# (c) James S Blachly, MD 2013
# 

# let print() be a function rather than statement
# ala python3
from __future__ import print_function

# python standard modules
import os
import sys
import time
import argparse #transitioned from getopt
import csv
import itertools
from collections import defaultdict
import gzip
import cPickle as pickle
import pdb   #pdb.set_trace()
import xml.etree.ElementTree as ET
import json

# nonstandard modules
import numpy as np
import pandas as pd
import HTSeq

# optional modules
try:
    import xlwt
except ImportError:
    print("Excel writer module xlwt not found; Microsoft Excel output disabled")

# mucor modules
import mucorfilters as mf
from variant import Variant
from mucorfeature import MucorFeature
import output
from config import Config
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
{0} [-h] | -g featurefile.gff -f feature_type [-u] -o <outputDir> <mutect001.txt mutect002.txt ... mutectNNN.txt>

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

def is_int(term):
    try:
        if int(term):
            return True
    except:
        return False 

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
    """"""
    
    config = Config()

    global filename2samples
    filename2samples = {}
    JD = json.load(open(json_config,'r'))
    config.featureType = JD['feature']
    config.outputDir = JD['outputDir']
    config.union = JD['union']
    config.fast = JD['fast']
    config.gff = JD['gff']
    config.outputFormats = JD['outputFormats']
    if str(JD['database']) == str("[u'[]']"):
        config.database = []
    else:
        config.database = JD['database']
    if str(JD['regions']):
        config.regions = JD['regions']
    else:
        config.regions = []
        
    config.filters = JD['filters']
    global database_switch
    database_switch = bool(config.database)
    global SnpEff_switch
    SnpEff_switch = False

    '''
    for i in JD['filters']:
        filters[i] = True      # Imagine filters as "ON/OFF", binary switches
    '''
    config.inputFiles = []
    config.samples = []
    for i in JD['samples']:
        config.samples.append(i['id'])
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
            config.inputFiles.append(j['path'])
    '''
    config.featureType = featureType
    config.outputDir = outputDir
    config.union = union
    config.fast = fast
    config.gff = gff
    config.database = database
    config.filters = filters
    config.inputFiles = inputFiles
    config.outputFormats = outputFormats
    '''

    return config #featureType, outputDir, union, fast, gff, database, filters, inputFiles, 

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

def parseRegionBed(regionfile, regionDictIn):
    regionDict = regionDictIn
    for line in open(str(regionfile),'r'):
        col = line.split("\t")
        chrom = col[0]
        start = col[1]
        end = col[2]
        regionDict[chrom].add((start,end))
    return regionDict

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

def parse_MuTectOUT(row, fieldId): 
    VF = row[fieldId['tumor_f']]
    DP = int(int(str(row[fieldId['t_ref_count']]).strip()) + int(str(row[fieldId['t_alt_count']]).strip()))
    position = int(row[fieldId['position']])
    '''
    MuTect_Annotations[tuple(( str(row[0]), position) )] = row[fieldId['dbsnp_site']]
    '''

    return VF, DP, position 
    
def parse_MuTectVCF(row, fieldId, header, fn): 
    j = 0
    for i in header:
        if str(i) not in ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'none']: # This line should detect if the sample id is in the line.
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
                MuTect_Annotations[tuple((str(row[0]), position))] = line.split("\t")[8]
                break
            else:
                continue
    '''

    return VF, DP, position 

def parse_SomaticIndelDetector(row, fieldId, header):
    j = 0
    '''
    #only works for samples with a dash (-) in them
    for i in header:
        if str('-') in str(i):
            tmpsampID = i
    '''
    # <solution to above> 
    #       assumes that sample ID is the final column in the header. always true? 
    #       if not always true, adopt the parse_mutect solution here as well
    tmpsampID = header[-1]
    # </solution to above>
    
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

def parse_HapCaller(row, fieldId, header):
    j = 0
    position = int(row[fieldId['POS']])
    for i in row[fieldId['FORMAT']].split(':'):
        if str(i) == "DP":
            DP = int(row[fieldId[header[-1]]].split(':')[j])
        if str(i) == "AD":
            AD = str(row[fieldId[header[-1]]].split(':')[j])
            if str(',') in AD:
                ref = int(AD.split(',')[0])
                alt = int(AD.split(',')[1])
                VF = float( float(alt)/(float(ref) + float(alt)) )
            else:
                abortWithMessage("Sample {0} may not have Haplotype Caller mutations with no ALT or VF".format(header[-1]))
        j += 1
    return VF, DP, position

def parse_FreeBayes(row, fieldId, header):
    j = 0
    position = int(row[fieldId['POS']])
    for i in row[fieldId['FORMAT']].split(':'):
        if str(i) == "DP":
            DP = int(row[fieldId[header[-1]]].split(':')[j])
        if str(i) == "RO":
            RO = int( str(row[fieldId[header[-1]]].split(':')[j]) )
        if str(i) == "AO":
            AO = int(  sum([ int(x) for x in str(row[fieldId[header[-1]]].split(':')[j]).split(',')]) )
        j += 1
    VF = float( float(AO)/float(AO + RO) )
    return VF, DP, position

def parse_MAF(row, fieldId, header):
    position = int(str(row[fieldId['Start_position']]).split('.')[0]) # case sensitive. what if, 'Start_Position' instead? 
    DP = int(str(row[fieldId['TTotCov']]).split('.')[0])
    VF = float( float(row[fieldId['TVarCov']])/float(DP) )
    chrom = str(row[fieldId['Chromosome']])
    ref =  str(row[fieldId['Reference_Allele']])
    alt = str(row[fieldId['Tumor_Seq_Allele2']])
    if ref == "-":
        ref = ""
    if alt == "-":
        alt = ""
    return VF, DP, position, chrom, ref, alt
'''
def parse_GVF(row, fieldId, header):
    position = 
'''
def filterRow(row, fieldId, filters, kind):
    '''
    returning True means this row will be filtered out [masked]
    returning False means the row will not be filtered out [pass filter]
    '''
    if str(kind) == "vcf":
        for rowFilter in str(row[fieldId['FILTER']]).split(';'):    ## VCF file format
            if rowFilter not in filters:
                return True
                break
    if str(kind) == "out":
        for rowFilter in str(row[fieldId['judgement']]).split(';'): ## MuTect '.out' file format
            if rowFilter not in filters:
                return True
                break
    return False
'''
def inRegion(chrom, start, end, region):
    region_chrom = region.split(':')[0]
    try:
        region_locs = region.split(':')[1]
        region_start = region_locs.split('-')[0]
        region_end = region_locs.split('-')[1]
        if str(chrom) == str(region_chrom) and int(start) >= int(region_start) and int(end) <= int(region_end):
            return True
    except IndexError:
        if str(chrom) == str(region_chrom):
            return True
    return False
'''

def inRegionDict(chrom, start, end, regionDict):
    if regionDict[chrom]:
        for locs in regionDict[chrom]:
            if locs[0] == 0 and locs[1] == 0:
                return True
            elif int(start) >= int(locs[0]) and int(end) <= int(locs[1]):
                return True
    return False

def skipThisIndel(ref, alt, knownFeatures, featureName, position, var):
    if len(var.ref) != len(var.alt):
        for kvar in knownFeatures[featureName].variants:
            if kvar.pos.pos == position and kvar.pos.chrom == var.pos.chrom and ( kvar.ref != var.ref or kvar.alt != var.alt ) and str(indelDelta(var.ref,var.alt)) == str(indelDelta(kvar.ref, kvar.alt)):
                return True, kvar.ref, kvar.alt
                break
    return False, '', ''

def indelDelta(ref, alt):
    ''' Detects the inserted or deleted bases of an indel'''
    if len(alt) > len(ref):             # insertion
        return alt.replace(ref,"",1)
    elif len(alt) < len(ref):           # deletion
        return ref.replace(alt,"",1)
    else:                               # SNP / MNP
        # there is already a SNP at this indel location 
        return ""

def parseVariantFiles(variantFiles, knownFeatures, gas, database, filters, regions): 
    # parse the variant files (VCF, muTect format)
    startTime = time.clock()

    # All variants stored in long (record) format
    # in a pandas dataframe
    varD = {'chr':[],'pos':[],'ref':[],'alt':[],'vf':[],'dp':[],'feature':[],'effect':[],'fc':[],'datab':[],'sample':[],'source':[]}
    
    if regions: # has the user specified any particular regions or region files to focus on?
        regionDict = defaultdict(set)
        for item in regions:
            if str(str(item).split('.')[-1]).lower() == 'bed':      # this item is a bed file
                regionDict = parseRegionBed(item, regionDict)
            elif str(str(item).split(':')[0]).startswith('chr'):    # this is a string 
                reg_chr = str(item.split(':')[0])
                try:
                    reg_str = str(str(item.split(':')[1]).split('-')[0])
                    reg_end = str(str(item.split(':')[1]).split('-')[1])
                except IndexError:                                  # represent whole chromosome regions [ex: chr2] by chrN:0-0 in the region dictionary   
                    reg_str = 0
                    reg_end = 0
                regionDict[reg_chr].add((reg_str, reg_end))

    print("\n=== Reading Variant Files ===")
    for fn in variantFiles:
        kind = str(str(fn).split('.')[-1].strip().lower())
        varFile = open(fn, 'rb')    # TO DO: error handling
        varReader = csv.reader(varFile, delimiter="\t")

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
        global HapCaller
        HapCaller = False
        global FreeBayes
        FreeBayes = False
        global MAF
        MAF = False
        global GVF
        GVF = False
        if str(fn.split('.')[-1].strip()).lower() == 'maf':
            MAF = True
            while str(row[0]).startswith('#'):
                row = varReader.next()
        if str(fn.split('.')[-1].strip()).lower() == 'gvf':
            GVF = True
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
            elif str('ID=HaplotypeCaller') in str(row):
                HapCaller = True
            elif str('freeBayes') in str(row):
                FreeBayes = True
            row = varReader.next()
        
        header = row
        if len(header) == 0: raise ValueError('Invalid header')
        fieldId = dict(zip(header, range(0, len(header))))

        # read coverage depth minimum cutoff; currently unusued
        #read_depth_min = 0
        # after reading the two header rows, read data
        for row in itertools.islice(varReader, None):
            if filterRow(row, fieldId, filters, kind):    # filter rows as they come in, to prevent them from entering the dataframe
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
                chrom = row[0]
                ref = row[3]
                alt = row[4]
            elif IonTorrent:
                VF, DP, position = parse_IonTorrent(row, fieldId, header)
                chrom = row[0]
                ref = row[3]
                alt = row[4]
            elif Mutect:
                VF, DP, position = parse_MuTectVCF(row, fieldId, header, fn) # , MuTect_Annotations)
                chrom = row[0]
                ref = row[3]
                alt = row[4]
            elif SomaticIndelDetector:
                VF, DP, position = parse_SomaticIndelDetector(row, fieldId, header)
                chrom = row[0]
                ref = row[3]
                alt = row[4]
            elif Mutector:
                VF, DP, position = parse_MuTectOUT(row, fieldId) # , MuTect_Annotations)
                chrom = row[0]
                ref = row[3]
                alt = row[4]
            elif Samtools:
                VF, DP, position = parse_SamTools(row, fieldId, header)
                chrom = row[0]
                ref = row[3]
                alt = row[4]
            elif VarScan:
                VF, DP, position = parse_VarScan(row, fieldId, header)
                chrom = row[0]
                ref = row[3]
                alt = row[4]
            elif HapCaller:
                VF, DP, position = parse_HapCaller(row, fieldId, header)
                chrom = row[0]
                ref = row[3]
                alt = row[4]
            elif FreeBayes:
                VF, DP, position = parse_FreeBayes(row, fieldId, header)
                chrom = row[0]
                ref = row[3]
                alt = row[4]
            elif MAF:
                VF, DP, position, chrom, ref, alt = parse_MAF(row, fieldId, header)

            else:
                abortWithMessage("{0} isn't a known data type:\nMiSeq, IonTorrent, SomaticIndelDetector, Samtools, VarScan, Haplotype Caller, or Mutect".format(fn))
            if is_int(chrom):
                chrom = str("chr" + str(chrom))
            if regions and not inRegionDict(chrom, int(position), int(position), regionDict ):
                continue
            var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=VF, dp=DP, eff=EFF.strip(';'), fc=FC.strip(';'))
            ###########################################
            # find bin for variant location
            resultSet = gas[ var.pos ]      # returns a set of zero to n IDs (e.g. gene symbols)
            if resultSet:                   # which I'll use as a key on the knownFeatures dict
                #print var.pos              # and each feature with matching ID gets allocated the variant
                #print(gas[ var.pos ])      # 
                for featureName in resultSet:
                    skipDupIndel, kvarref, kvaralt = skipThisIndel(var.ref, var.alt, knownFeatures, featureName, position, var)
                    if bool(skipDupIndel):
                        # Sanity check to see what indels are being overwritten by existing vars
                        #print(var.pos.chrom, var.pos.pos, var.ref, var.alt, kvarref, kvaralt)
                        var.ref = kvarref
                        var.alt = kvaralt
                    
                    knownFeatures[featureName].variants.add(var)
            
            # Descriptive variable names
            chr = chrom
            pos = int(position)
            #ref = row[3]
            #alt = row[4]
            vf = VF     # to do, change (karl capitalizes)
            dp = DP     # to do, change (karl capitalizes)
            feature = ', '.join( gas[ var.pos ] )   # join with comma to handle overlapping features
            effect = EFF
            fc = FC
            #sample = fn.split('/')[-1]      # to do, will need to come from JSON config
            sample = filename2samples[str(fn.split('/')[-1])]
            source = fn.split('/')[-1]
            rowAnnotated, rowAnnotation = isAnnotatedSNP(var, database)
            datab = str(rowAnnotation)

            '''
            if ( str('NM_001290815') in str(feature) or str(feature) == "NM_001290815" ) and str(chr) == "chr19" and int(pos) == 57110490:
                pdb.set_trace()
                # for i in knownFeatures['NM_001290815'].variants: print(i.pos.pos)
            '''


            '''
            if Mutect or Mutector:
                datab = MuTect_Annotations[(chr, pos)]
            '''
            '''
            if skipThisIndel(ref, alt):
                for var in knownFeatures[feature]
            '''
            # build dict to insert
            #columns=('chr','pos','ref','alt','vf','dp','gene','effect','sample','source')
            vardata = dict(zip( ['chr','pos','ref','alt','vf','dp','feature','effect','fc','datab','sample','source'], \
                                [ chr , pos , ref , alt , vf , dp , feature , effect , fc , datab , sample , source ])) 
            # add to variants data frame dictionary
            for key in vardata.keys():
                varD[key].append(vardata[key])

        totalTime = time.clock() - startTime
        print("{0:02d}:{1:02d}\t{2}".format(int(totalTime/60), int(totalTime % 60), fn))
    # Transform data frame dictionary into pandas DF. Major speed increase relative to appending the DF once per variant
    # Removing columns from the following 'columns' list will mask them from output in allvars.txt
    varDF = pd.DataFrame(varD, columns=['chr','pos','ref','alt','vf','dp','feature','effect','fc','datab','sample','source'])
    # Clean up variant dataframe a little
    # position should be integer, not float
    varDF.pos = varDF.pos.astype(int)
    return varDF, knownFeatures, gas 

def printRunInfo(config, outputDirName):
    try: 
        ofRunInfo = open(outputDirName + "/run_info.txt", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))
    # =========================
    # run_info.txt
    #
    ofRunInfo.write(Info.versionInfo + "\n")
    ofRunInfo.write("{0}\n\n".format(time.ctime() ) )
    ofRunInfo.write(str(config))
    
    '''
    ofRunInfo.write("Variants Pre-filter: \n")
    ofRunInfo.write("        Post-filter: \n")
    '''
    ofRunInfo.close()

def printCounts(outputDirName, knownFeatures):
    try:
        ofCounts = open(outputDirName + "/counts.txt", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))

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
            ofCounts.write(feature.name + "\t")

            ofCounts.write(str(len(knownFeatures[feature.name].variants)) + "\t")
            
            ofCounts.write(str(knownFeatures[feature.name].weightedVariants()) + "\t")
            
            ft = knownFeatures[feature.name]
            avgWt = float(ft.weightedVariants() / float(ft.numVariants()) )
            ofCounts.write(str(avgWt) + "\t")

            ofCounts.write(str(knownFeatures[feature.name].numUniqueVariants()) + "\t")
            
            ofCounts.write(str(knownFeatures[feature.name].numUniqueSamples()) + "\n")

            nrow += 1
    
    print("\t{0}: {1} rows".format(ofCounts.name, nrow))
    ofCounts.close()

def printVariantDetails(outputDirName, knownFeatures, varDF, total):
    try:
        ofVariantDetails = open(outputDirName + "/variant_details.txt", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))

    # =========================================================
    # variant_details.txt
    #

    ofVariantDetails.write('Feature\tContig\tPos\tRef\tAlt\tVF\tDP\t') ######## Karl Modified ##############
    if SnpEff_switch:
        ofVariantDetails.write('Effect\tFC\t')
    if database_switch:
        ofVariantDetails.write('Annotation\t')
    ofVariantDetails.write('Count\tFrequency\tSource\n')

    masterList = list(knownFeatures.values())
    sortedList = sorted(masterList, key=lambda k: k.numVariants(), reverse=True)
    nrow = 0

    for feature in sortedList:
        if knownFeatures[feature.name].variants:
            for var in knownFeatures[feature.name].uniqueVariants():
                
                #if not isAnnotatedSNP(snps, tuple((var.pos.chrom,var.pos.pos))):  ##### KARL ADDED ######
                ofVariantDetails.write(feature.name + "\t")
                ofVariantDetails.write(var.pos.chrom + "\t")
                ofVariantDetails.write(str(var.pos.pos) + "\t")
                ofVariantDetails.write(var.ref + "\t")
                ofVariantDetails.write(var.alt + "\t")
                ofVariantDetails.write(var.frac + "\t")
                ofVariantDetails.write(var.dp + "\t")
                
                if SnpEff_switch:
                    if str(var.eff) != str(''):
                        ofVariantDetails.write(str([ x for x in set(str(var.eff).split(', '))]).replace("''","").replace(", ","").strip(']').strip('[').strip("'") + "\t")
                    else:
                        ofVariantDetails.write(str('?') + "\t")
                    if str(var.fc) != str(''):
                        ofVariantDetails.write(str([ x for x in set(str(var.fc).split(', '))]).replace("''","").replace(", ","").strip(']').strip('[').strip("'") + "\t")
                    else:
                        ofVariantDetails.write(str('?') + "\t")
                if database_switch:
                    if len(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()) == 1:
                        ofVariantDetails.write(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()[0] + "\t") 
                    else:
                        tempdb = ""
                        for i in [x.strip(']').strip('[') for x in str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()).replace(',','').split(' ')]:
                            if str(i) not in str(tempdb):
                                tempdb += str(i) + ","
                        ofVariantDetails.write(tempdb.strip(',') + "\t")
                ofVariantDetails.write(str(len(var.source.split(','))) + "\t")
                ofVariantDetails.write(str(float(len(var.source.split(',')))/float(total)) + "\t" )
                ofVariantDetails.write(var.source + "\n")
                nrow += 1
                
        ########### Karl added bed file output here #############
    ofVariantDetails.close()
    print("\t{0}: {1} rows".format(ofVariantDetails.name, nrow))

def printVariantDetailsXLS(outputDirName, knownFeatures, varDF, total):
    '''
    try:
        ofVariantDetails = open(outputDirName + "/variant_details.txt", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))
    '''
    # =========================================================
    # variant_details.txt
    #
    workbook = xlwt.Workbook()
    ofVariantDetails = workbook.add_sheet('Variant Details')
    nrow = 0
    ncol = 0
    for column in ['Feature', 'Contig', 'Pos', 'Ref', 'Alt', 'VF', 'DP']:
        ofVariantDetails.write(nrow, ncol, column)
        ncol += 1
    if SnpEff_switch:
        for column in ['Effect', 'FC']:
            ofVariantDetails.write(nrow, ncol, column)
            ncol += 1
    if database_switch:
        ofVariantDetails.write(nrow, ncol, 'Annotation')
        ncol += 1
    for column in ['Count','Frequency','Source']:
        ofVariantDetails.write(nrow, ncol, column)
        ncol += 1
    nrow += 1
    masterList = list(knownFeatures.values())
    sortedList = sorted(masterList, key=lambda k: k.numVariants(), reverse=True)
    

    for feature in sortedList:
        if knownFeatures[feature.name].variants:
            for var in knownFeatures[feature.name].uniqueVariants():
                ncol = 0
                ofVariantDetails.write(nrow, ncol, str(feature.name))
                ncol += 1
                ofVariantDetails.write(nrow, ncol, str(var.pos.chrom))
                ncol += 1
                ofVariantDetails.write(nrow, ncol, int(var.pos.pos))
                ncol += 1
                ofVariantDetails.write(nrow, ncol, str(var.ref))
                ncol += 1
                ofVariantDetails.write(nrow, ncol, str(var.alt))
                ncol += 1
                ofVariantDetails.write(nrow, ncol, str(var.frac))
                ncol += 1
                ofVariantDetails.write(nrow, ncol, str(var.dp))
                ncol += 1
                
                if SnpEff_switch:
                    if str(var.eff) != str(''):
                        ofVariantDetails.write(nrow, ncol, str(str([ x for x in set(str(var.eff).split(', '))]).replace("''","").replace(", ","").strip(']').strip('[').strip("'")))
                        ncol += 1
                    else:
                        ofVariantDetails.write(nrow, ncol, str('?'))
                        ncol += 1
                    if str(var.fc) != str(''):
                        ofVariantDetails.write(nrow, ncol, str(str([ x for x in set(str(var.fc).split(', '))]).replace("''","").replace(", ","").strip(']').strip('[').strip("'")))
                        ncol += 1
                    else:
                        ofVariantDetails.write(nrow, ncol, str('?'))
                        ncol += 1
                if database_switch:
                    if len(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()) == 1:
                        ofVariantDetails.write(nrow, ncol, str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()[0]) )
                        ncol += 1
                    else:
                        tempdb = ""
                        for i in [x.strip(']').strip('[') for x in str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()).replace(',','').split(' ')]:
                            if str(i) not in str(tempdb):
                                tempdb += str(i) + ","
                        ofVariantDetails.write(nrow, ncol, str(tempdb.strip(',')))
                        ncol += 1
                ofVariantDetails.write(nrow, ncol, int(str(len(var.source.split(',')))))
                ncol += 1
                ofVariantDetails.write(nrow, ncol, float(float(len(var.source.split(',')))/float(total))  )
                ncol += 1
                ofVariantDetails.write(nrow, ncol, str(var.source))
                ncol += 1
                nrow += 1
                
        ########### Karl added bed file output here #############
    workbook.save(outputDirName + '/variant_details.xls')
    print("\t{0}: {1} rows".format(str(outputDirName + '/variant_details.xls'), nrow))

def printLongVariantDetailsXLS(outputDirName, knownFeatures, varDF, total):
    '''
    try:
        ofVariantDetails = open(outputDirName + "/variant_details.txt", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))
    '''
    # =========================================================
    # variant_details.txt
    #
    workbook = xlwt.Workbook()
    ofLongVariantDetails = workbook.add_sheet('Long Variant Details')
    nrow = 0
    ncol = 0
    for column in ['Feature', 'Contig', 'Pos', 'Ref', 'Alt', 'VF', 'DP']:
        ofLongVariantDetails.write(nrow, ncol, column)
        ncol += 1
    if SnpEff_switch:
        for column in ['Effect', 'FC']:
            ofLongVariantDetails.write(nrow, ncol, column)
            ncol += 1
    if database_switch:
        ofLongVariantDetails.write(nrow, ncol, 'Annotation')
        ncol += 1
    for column in ['Count','Frequency','Source']:
        ofLongVariantDetails.write(nrow, ncol, column)
        ncol += 1
    nrow += 1
    masterList = list(knownFeatures.values())
    sortedList = sorted(masterList, key=lambda k: k.numVariants(), reverse=True)
    

    for feature in sortedList:
        if knownFeatures[feature.name].variants:
            for var in knownFeatures[feature.name].uniqueVariants():
                #pdb.set_trace()
                for j in range(len(var.dp.split(','))):
                    ncol = 0
                    ofLongVariantDetails.write(nrow, ncol, str(feature.name))
                    ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, str(var.pos.chrom))
                    ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, int(var.pos.pos))
                    ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, str(var.ref))
                    ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, str(var.alt))
                    ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, float((var.frac).split(', ')[j]))
                    ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, int((var.dp).split(', ')[j]))
                    ncol += 1
                    
                    if SnpEff_switch:
                        if str(var.eff) != str(''):
                            ofLongVariantDetails.write(nrow, ncol, str(str([ x for x in set(str(var.eff).split(', '))]).replace("''","").replace(", ","").strip(']').strip('[').strip("'")))
                            ncol += 1
                        else:
                            ofLongVariantDetails.write(nrow, ncol, str('?'))
                            ncol += 1
                        if str(var.fc) != str(''):
                            ofLongVariantDetails.write(nrow, ncol, str(str([ x for x in set(str(var.fc).split(', '))]).replace("''","").replace(", ","").strip(']').strip('[').strip("'")))
                            ncol += 1
                        else:
                            ofLongVariantDetails.write(nrow, ncol, str('?'))
                            ncol += 1
                    if database_switch:
                        if len(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()) == 1:
                            ofLongVariantDetails.write(nrow, ncol, str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()[0]) )
                            ncol += 1
                        else:
                            tempdb = ""
                            for i in [x.strip(']').strip('[') for x in str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()).replace(',','').split(' ')]:
                                if str(i) not in str(tempdb):
                                    tempdb += str(i) + ","
                            ofLongVariantDetails.write(nrow, ncol, str(tempdb.strip(',')))
                            ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, int(len(var.source.split(','))))
                    ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, float(float(len(var.source.split(',')))/float(total))  )
                    ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, str(var.source.split(', ')[j]))
                    ncol += 1
                    nrow += 1
                
        ########### Karl added bed file output here #############
    workbook.save(outputDirName + '/long_variant_details.xls')
    print("\t{0}: {1} rows".format(str(outputDirName + '/long_variant_details.xls'), nrow))

def printVariantBed(outputDirName, knownFeatures):
    try:
        ofVariantBeds = open(outputDirName + "/variant_locations.bed", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))
    # =========================================================
    # variant_details.txt
    #
    masterList = list(knownFeatures.values())
    sortedList = sorted(masterList, key=lambda k: k.numVariants(), reverse=True)
    nrow = 0

    for feature in sortedList:
        if knownFeatures[feature.name].variants:
            for var in knownFeatures[feature.name].uniqueVariants():
                ofVariantBeds.write(var.pos.chrom + "\t")
                ofVariantBeds.write(str(var.pos.pos - 1) + "\t")
                ofVariantBeds.write(str(var.pos.pos) + "\t")
                #ofVariantBeds.write(str(snps[tuple((var.pos.chrom,var.pos.pos))]) + "\t") # used to write dbSNP name in bed name field
                ofVariantBeds.write("\n")
                nrow += 1
    ofVariantBeds.close()
    print("\t{0}: {1} rows".format(ofVariantBeds.name, nrow))

def getMetricsVCF(var, sourcefile):
    n = 0
    dp = 0
    vf = 0
    for source in str(var.source).split(', '):
        if source == sourcefile:
            dp = str(var.dp).split(', ')[n]
            vf = str(var.frac).split(', ')[n]
        else:
            n += 1
    if dp and vf:
        return dp, vf
    elif not dp and not vf:
        return 0, 0
    else:
        abortWithMessage("There should be a depth AND variant frequency for this mutation. One is missing.")


def printBigVCF(outputDirName, knownFeatures, varDF, inputFiles):
    try:
        ofBigVCF = open(outputDirName + "/variant_calls.vcf", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))
    # =========================================================
    # variant_calls.vcf
    #    
    ofBigVCF.write("##fileformat=VCFv4.1\n")
    ofBigVCF.write('##FORMAT=<ID=DP,Type=Float,Description="Depth of read coverage">\n')
    ofBigVCF.write('##FORMAT=<ID=VF,Type=Float,Description="Fraction of reads displaying alternative allele">\n')
    ofBigVCF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
    for inFile in inputFiles:
        ofBigVCF.write(inFile.split('/')[-1] + "\t")
    ofBigVCF.write("\n")
    masterList = list(knownFeatures.values())
    sortedList = sorted(masterList, key=lambda k: k.numVariants(), reverse=True)
    nrow = 0

    for feature in sortedList:
        if knownFeatures[feature.name].variants:
            for var in knownFeatures[feature.name].uniqueVariants():
                ofBigVCF.write(str(var.pos.chrom) + "\t" + str(var.pos.pos) + "\t")
                if database_switch:
                    if len(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()) == 1:
                        ofVariantDetails.write(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()[0] + "\t") 
                    else:
                        tempdb = ""
                        for i in [x.strip(']').strip('[') for x in str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom))].datab.unique()).replace(',','').split(' ')]:
                            if str(i) not in str(tempdb):
                                tempdb += str(i) + ";"
                        ofVariantDetails.write(tempdb.strip(';') + "\t")
                else:
                    ofBigVCF.write(str(".") + "\t")
                if var.ref:
                    ofBigVCF.write(str(var.ref) + "\t") 
                else:
                    ofBigVCF.write("." + "\t")
                if var.alt:
                    ofBigVCF.write(str(var.alt) + "\t")
                else:
                    ofBigVCF.write("." + "\t")

                ofBigVCF.write("." + "\t" + "." + "\t" + "." + "\t" + "DP:VF" + "\t")
                for source in inputFiles:
                    dp, vf = getMetricsVCF(var, source.split('/')[-1])
                    ofBigVCF.write(str(dp) + ":" + str(vf) + "\t")
                ofBigVCF.write("\n")
                nrow += 1
    ofBigVCF.close()
    print("\t{0}: {1} rows".format(ofBigVCF.name, nrow))


def printOutput(config, outputDirName, outputFormats, knownFeatures, gas, varDF): ######## Karl Modified ##############
    '''Output statistics and variant details to the specified output directory.'''

    startTime = time.clock()
    print("\n=== Writing output files to {0}/ ===".format(outputDirName))

    # TODO: pull output options from json file (config class) and put if statements here to selectively output
    outputFormatsDict = defaultdict(bool)
    for format in outputFormats:
        outputFormatsDict[format] = bool(True)
    total = len(set(config.samples))                # used in frequency calculation. config.samples will use sample count as the denominator. samples.inputFiles will use file count
    printRunInfo(config, outputDirName)
    if outputFormatsDict['counts']:
        printCounts(outputDirName, knownFeatures)
    if outputFormatsDict['txt']:
        printVariantDetails(outputDirName, knownFeatures, varDF, total)
    if outputFormatsDict['bed']:
        printVariantBed(outputDirName, knownFeatures)
    if outputFormatsDict['xls'] and 'xlwt' in sys.modules:
        printVariantDetailsXLS(outputDirName, knownFeatures, varDF, total)
    if outputFormatsDict['default']:
        if 'xlwt' in sys.modules: printVariantDetailsXLS(outputDirName, knownFeatures, varDF, total)
        printVariantBed(outputDirName, knownFeatures)
        printCounts(outputDirName, knownFeatures)
    if outputFormatsDict['long'] and 'xlwt' in sys.modules:
        printLongVariantDetailsXLS(outputDirName, knownFeatures, varDF, total)
    if outputFormatsDict['vcf']:
        printBigVCF(outputDirName, knownFeatures, varDF, config.inputFiles)

    totalTime = time.clock() - startTime
    print("\tTime to write: {0:02d}:{1:02d}".format(int(totalTime/60), int(totalTime % 60)))
    return 0

def main():
    print(Info.logo)
    print(Info.versionInfo)

    print("\n=== Run Info ===")
    print("\t{0}".format(time.ctime() ) )
    print()

    config = parseJSON(sys.argv[1])

    if not os.path.exists(config.gff):
        abortWithMessage("Could not find GFF file {0}".format(config.gff))
    if os.path.exists(config.outputDir) and os.listdir(config.outputDir):
        abortWithMessage("The directory {0} already exists and contains output. Will not overwrite.".format(config.outputDir))
    elif not os.path.exists(config.outputDir):
        try:
            os.makedirs(config.outputDir)
        except:
            abortWithMessage("Error when creating output directory {0}".format(config.outputDir))

    # check that all specified variant files exist
    for fn in config.inputFiles:
        if not os.path.exists(fn):
            abortWithMessage("Could not find variant file: {0}".format(fn))

    knownFeatures, gas = parseGffFile(str(config.gff), str(config.featureType), bool(config.fast))

    varDF, knownFeatures, gas = parseVariantFiles(list(config.inputFiles), knownFeatures, gas, config.database, config.filters, config.regions)
    printOutput(config, str(config.outputDir), config.outputFormats, knownFeatures, gas, varDF) ######## Karl Modified ##############
    
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
    ##varDF.sort(columns=['feature','pos'], inplace=True)
    ##varDF.replace('', np.nan, inplace=True)
    ##varDF.to_csv(outputDir + '/allvars.txt', sep="\t", na_rep='?', index=False)
    
    # FUTURE: all output taken care of with below
    # ow = output.writer()
    # for format in config.outputFormats
    #   ow.write(df, format, config.outputDir) #<-outputDir or whtever, double check

    # >>>>>>>>> JAMES' OUTPUT FUNCTION  <<<<<<<<<<<
    # ow = output.writer()
    # ow.write(varDF, "long", outputDir)

    # pretty print newline before exit
    print()


if __name__ == "__main__":
    if sys.hexversion < 0x02070000:
        raise RuntimeWarning("mucor should be run on python 2.7.0 or greater.")
    main()
