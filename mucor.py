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
import argparse
import csv
import itertools
from collections import defaultdict
import gzip
import cPickle as pickle
import pdb                             # Needed for pdb.set_trace()
import xml.etree.ElementTree as ET
import json

# nonstandard modules
import numpy as np
import pandas as pd
from pandas import ExcelWriter
import HTSeq

# optional modules
try:
    import xlwt
except ImportError:
    print("Excel writer module xlwt not found; Microsoft Excel output disabled")
try:
    import tabix
except ImportError:
    print("Tabix module not found; database features disabled")

# mucor modules
import mucorfilters as mf
from variant import Variant
from mucorfeature import MucorFeature
from inputs import *
import output
from config import Config
from databases import dbLookup 

class Info:
    '''Program info: logo, version, and usage'''
    logo = """
 __    __     __  __     ______     ______     ______    
/\ "-./  \   /\ \/\ \   /\  ___\   /\  __ \   /\  == \   
\ \ \-./\ \  \ \ \_\ \  \ \ \____  \ \ \/\ \  \ \  __<   
 \ \_\ \ \_\  \ \_____\  \ \_____\  \ \_____\  \ \_\ \_\ 
  \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_/ /_/ 
                                  
"""
    version = "0.9"
    versionInfo = "mucor version {0}\nJames S Blachly, MD\nKarl W Kroll, BS".format(version)
    # usage needs to be updated, or eliminated if argparse can replace this function. 
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

def throwWarning(message, help = False):
    print("*** WARNING: " + message + " ***")
    return

def constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures):
    gas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for feature in itertools.islice(gffFile, 0, None):
        # Nonstandard contigs (eg chr17_ctg5_hap1, chr19_gl000209_random, chrUn_...)
        # must be specifically excluded, otherwise you will end up with exception
        # ValueError: start is larger than end
        # due to duplicated gene symbols
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
    """
    Import the JSON config file from mucor_config.py. 
    Reads the config file into a dictionary, then writes each dictionary entry into the respective Config class position.
    """
    config = Config()

    global filename2samples
    filename2samples = {}
    # load the json file into a dictionary, JD
    JD = json.load(open(json_config,'r'))
    # write dictionary values into a more human-friendly config class
    config.featureType = JD['feature']
    config.outputDir = JD['outputDir']
    config.union = JD['union']
    config.fast = JD['fast']
    config.gff = JD['gff']
    config.outputFormats = JD['outputFormats']
    if JD['databases']:
        config.databases = JD['databases']
    if str(JD['regions']):
        config.regions = JD['regions']
    else:
        config.regions = []
    # comma separated list of acceptable VCF filter column values
    config.filters = JD['filters']

    global database_switch
    # if any databases are defined, switch the database_switch ON (True)
    # if no databases are defined, config.database will be [] and database_switch will be OFF (False)
    if 'tabix' in sys.modules:
        database_switch = bool(config.databases)
    else:
        database_switch = bool(False)
    
    global SnpEff_switch
    # works similar to database switch above.
    SnpEff_switch = False

    config.inputFiles = []
    config.samples = []
    for i in JD['samples']:
        config.samples.append(i['id'])
        for j in i['files']:
            filename = str(j['path']).split('/')[-1]
            filename2samples[filename] = i['id']
            if not SnpEff_switch and str(j['type']) == str('vcf') and bool(j['snpeff']) == bool(True):
                SnpEff_switch = bool(True)
            config.inputFiles.append(j['path'])

    return config 

def parseGffFile(gffFileName, featureType, fast):
    '''
    Parse the GFF/GTF file. Return tuple (knownFeatures, GenomicArrayOfSets)
    Haplotype contigs are explicitly excluded because of a coordinate crash (begin > end)
    '''
    
    # TO DO: command line flag should indicate that variants in INTRONS are counted
    # This is called --union, see below
    
    startTime = time.clock()
    print("\n=== Reading GFF/GTF file {0} ===".format(gffFileName))
    
    scripts = "/".join(os.path.realpath(__file__).split('/')[:-1])
    annotFileName = ".".join(gffFileName.split('/')[-1].split('.')[:-1])
    # pickled file is created for specific combinations of gff annotation and feature type. See below for more details. **
    archiveFilePath = str("/") + str(scripts).strip('/') + str("/") + str(annotFileName) + str('_') + str(featureType) + str('.p')
    print(gffFileName)
    gffFile = HTSeq.GFF_Reader(gffFileName)
    
    knownFeatures = {}

    duplicateFeatures = set()

    # gas - GenomicArrayOfSets.
    # typecode always 'O' (object) for GenomicArrayOfSets
    # UNstranded -- VCF and muTect output always report on + strand,
    # but the GenomicArray must be unstranded because the GFF /is/ strand-specific,
    # and if I manually coded all GenomicIntervals read from the VCF or muTect file as '+',
    # then no genes on the - strand would have variants binned to them

    # Fast boolean determines whether the user wants to load a pickle file annotation that they made in a previous run
    # The pickle file will load the genomic array of sets, known features, and duplicate features. 
    # ** These three items will change depending on the supplied gff annotation AND the feature selected. 
    #    Thus, each pickle file will be for a specific combination of gff and feature, hence the naming convention.
    if fast:
        # user wants to use pickled annotations
        if os.path.exists( archiveFilePath ):
            # using the existing, pickled annotation
            print("Opening annotation archive: " + str(archiveFilePath))
            pAnnot = open(archiveFilePath, 'rb')
            gas = pickle.load(pAnnot)
            knownFeatures = pickle.load(pAnnot)
            duplicateFeatures = pickle.load(pAnnot)
            pAnnot.close()
        else:
            # no pickled annotation exists for this combination of gff and feature; creating it in the same directory as this mucor.py file
            print("Cannot locate annotation archive for " + str(gffFileName.split('/')[-1]) + str(" w/ ") + str(featureType))
            print("   Reading in annotation and saving archive for faster future runs") 
            gas, knownFeatures, duplicateFeatures = constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures)
            archiveOut = open(archiveFilePath, 'wb')
            pickle.dump(gas, archiveOut, -1) ### pickle feature only works with full annotation files
            pickle.dump(knownFeatures, archiveOut, -1)
            pickle.dump(duplicateFeatures, archiveOut, -1)
            archiveOut.close()
    if not fast:
    # ignore pickles function entirely. Won't check for it and won't attempt to create it
        gas, knownFeatures, duplicateFeatures = constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures)

    if duplicateFeatures:
        print("*** WARNING: {0} {1}s found on more than one contig".format(len(duplicateFeatures), featureType))

    totalTime = time.clock() - startTime
    print("{0} sec\t{1} found:\t{2}".format(int(totalTime), featureType, len(knownFeatures)))

    return knownFeatures, gas

def parseRegionBed(regionfile, regionDictIn):
    ''' 
    Read through the supplied bed file and add the rows to the input region dictionary before returning it. Appends to the input dict, rather than overwriting it.
    '''
    regionDict = regionDictIn
    for line in open(str(regionfile),'r'):
        col = line.split("\t")
        chrom = col[0]
        start = col[1]
        end = col[2]
        regionDict[chrom].add((start,end))
    return regionDict

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

def inRegionDict(chrom, start, end, regionDict):
    '''Checks the given mutation location to see if it is in the dictionary of regions'''
    if regionDict[chrom]: # are there any regions of interest in the same chromosome as this mutation?
        for locs in regionDict[chrom]: # breaking the list of regions according to chromosome should significantly decrease the number of comparisons necessary 
            if locs[0] == 0 and locs[1] == 0: # chrN:0-0 is used to define an entire chromosome as a region of interest. 
                return True
            elif int(start) >= int(locs[0]) and int(end) <= int(locs[1]):
                return True
    return False

def skipThisIndel(var, knownFeatures, featureName):
    '''
    Input a mutation and the current list of known mutations and features
    If this mutation already exists in another form, returns True, along with the existing ref and alt.
        Existing in another form is defined as 2 mutations that occur in the same position and have different ref/alt, but are functionally equivalent
        Ex: ref/alt: A/AT compared to ref/alt: ATT/ATTT
        The refs and the alts are different, but both mutations represent the same T insertion
    '''
    if len(var.ref) != len(var.alt):
        for kvar in knownFeatures[featureName].variants: # kvar --> known variant
            if kvar.pos.pos == var.pos.pos and kvar.pos.chrom == var.pos.chrom and ( kvar.ref != var.ref or kvar.alt != var.alt ) and str(indelDelta(var.ref,var.alt)[0]) == str(indelDelta(kvar.ref, kvar.alt)[0]) and str(indelDelta(var.ref,var.alt)[1]) == str(indelDelta(kvar.ref, kvar.alt)[1]):
                return (kvar.ref, kvar.alt)
                break
    return tuple()

def indelDelta(ref, alt):
    ''' Detects the inserted or deleted bases of an indel'''
    if len(alt) > len(ref):             # insertion
        return [ alt.replace(ref,"",1), str("INS") ]
    elif len(alt) < len(ref):           # deletion
        return [ ref.replace(alt,"",1), str("DEL") ]
    else:                               # SNP / MNP
        # there is already a SNP at this indel location 
        return [ "", "" ]

def parseVariantFiles(variantFiles, knownFeatures, gas, databases, filters, regions): 
    '''
    Read in all input files
    Log mutations in the variant dataframe 
    '''
    startTime = time.clock()

    # All variants stored in long (record) format
    # in a pandas dataframe
    #
    # However, it is MUCH faster to initially store them in a python Dict
    # Then convert to the pandas DF at the end
    varD = defaultdict(list) # {'chr':[],'pos':[],'ref':[],'alt':[],'vf':[],'dp':[],'feature':[],'effect':[],'fc':[],'datab':[],'sample':[],'source':[]}
    
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
        
        try:
            varFile = open(fn, 'rb')
        except IOError:
            # file does not exist or cannot be opened
            throwWarning('{} could not be opened'.format(fn))
            continue    # next fn in variantFiles

        varReader = csv.reader(varFile, delimiter="\t")

        try:
            row = varReader.next()
        except StopIteration:
            # a file was empty (i.e. no first row to read)
            print("Empty file {}".format(fn))
            continue    # next fn in variantFiles

        # declare all data types False, then reassign to True as they are encountered.
        # only one should be true at a time for any given file (fn)
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

        # check filename for MAF or GVF
        if str(fn.split('.')[-1].strip()).lower() == 'maf':
            MAF = True
            while str(row[0]).startswith('#'):
                row = varReader.next()
        if str(fn.split('.')[-1].strip()).lower() == 'gvf':
            GVF = True

        # check VCF file headers for the variant calling software used
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

        # after reading the two header rows, read data
        for row in itertools.islice(varReader, None):
            if filterRow(row, fieldId, filters, kind):  # filter rows as they come in, to prevent them from entering the dataframe
                continue                                # this allows us to print the dataframe directly and have consistent output with variant_details.txt, etc.

            effect = ""
            fc = ""
            muts = []
            loca = []
            try:
                for eff in row[fieldId['INFO']].split(';'):
                    if eff.startswith('EFF='):
                        for j in eff.split(','):
                            muts.append(str(j.split('|')[3]))
                            loca.append(str(j.split('(')[0]).replace('EFF=',''))
                for mut in set(muts):
                    if str(mut) != "":
                        effect += str(mut) + ";"
                for loc in set(loca):
                    if str(loc) != "":
                        fc += str(loc) + ";"
            except KeyError:
                pass

            # prime the InputParser object for this row
            InputParser.row = row
            InputParser.fieldId = fieldId
            InputParser.header = header
            InputParser.fn = fn
            InputParser.eff = effect
            InputParser.fc = fc 

            # parse the row of data, depending on the type of data
            if MiSeq:
                var = parse_MiSeq(InputParser)                
            elif IonTorrent:
                var = parse_IonTorrent(InputParser)                
            elif Mutect:
                var = parse_MuTectVCF(InputParser)                
            elif SomaticIndelDetector:
                var = parse_SomaticIndelDetector(InputParser)
            elif Mutector:
                var = parse_MuTectOUT(InputParser)
            elif Samtools:
                var = parse_SamTools(InputParser)
            elif VarScan:
                var = parse_VarScan(InputParser)
            elif HapCaller:
                var = parse_HapCaller(InputParser)
            elif FreeBayes:
                var = parse_FreeBayes(InputParser)
            elif MAF:
                var = parse_MAF(InputParser)
            else:
                abortWithMessage("{0} isn't a known data type:\nMiSeq, IonTorrent, SomaticIndelDetector, Samtools, VarScan, Haplotype Caller, or Mutect".format(fn))
            if regions and not inRegionDict(chrom, int(position), int(position), regionDict ):
                continue

            ###########################################
            # find bin for variant location
            resultSet = gas[ var.pos ]      # returns a set of zero to n IDs (e.g. gene symbols)
            if resultSet:                   # which I'll use as a key on the knownFeatures dict
                #                           # and each feature with matching ID gets allocated the variant
                for featureName in resultSet:
                    kvar = skipThisIndel(var, knownFeatures, featureName)
                    if bool(kvar):
                        # Sanity check to see what indels are being overwritten by existing vars
                        var.ref = kvar[0]
                        var.alt = kvar[1]
                    
                    knownFeatures[featureName].variants.add(var)
            
            # Descriptive variable names
            chr = var.pos.chrom
            pos = int(var.pos.pos)
            ref = var.ref
            alt = var.alt 
            vf = var.frac
            dp = var.dp
            features = ', '.join( gas[ var.pos ] )   # join with comma to handle overlapping features
            effect = var.eff
            fc = var.fc
            sample = filename2samples[str(fn.split('/')[-1])]
            source = os.path.split(fn)[1]
            dbEntries = dbLookup(var, databases)
            count = 0 # initialize this database column now to save time later

            # build dict to insert
            # Step 1: define the columns and their values
            # Step 1b.The database columns are variable, ranging from zero to many
            #         and must be inserted dynamically
            # Step 2. zip the column names and column values together into a dictionary
            # Step 3. Add this round to the master variants data frame dictinoary

            for feature in features.split(', '):
                columns = ['chr','pos','ref','alt','vf','dp','feature','effect','fc']
                values  = [ chr,  pos,  ref,  alt,  vf,  dp,  feature,  effect,  fc]

                for dbName in sorted(dbEntries.keys()):
                    columns.append(dbName)
                    values.append(dbEntries[dbName])
                columns += ['count', 'sample','source']
                values  += [ count,   sample,  source ]

                vardata = dict(zip( columns, values ))
                for key in vardata.keys():
                    varD[key].append(vardata[key])

        totalTime = time.clock() - startTime
        print("{0:02d}:{1:02d}\t{2}".format(int(totalTime/60), int(totalTime % 60), fn))
    
    # Transform data frame dictionary into pandas DF. Major speed increase relative to appending the DF once per variant
    # Removing columns from the following 'columns' list will mask them from output
    varDF = pd.DataFrame(varD, columns=columns)
    print("Finished reading\nAssessing total variant counts")
    varDF['count'] = varDF.groupby(['chr', 'pos', 'ref', 'alt', 'feature'])['sample'].transform(len)
    # Clean up variant dataframe a little
    # position should be integer, not float
    varDF.pos   = varDF.pos.astype(int)
    
    return varDF, knownFeatures, gas 

def getCountsAndFrequency(varDF, total):
    mutDF = varDF.drop_duplicates(subset=['chr','pos','alt','feature'])[['chr','pos','alt','feature']] # makes a new dataframe of chr, pos, and alt; all you need to identify unique mutations
    for i in dedupvart.iterrows(): 
        count = len(varDF[(varDF.chr == i[1][0]) & (varDF.pos == i[1][1]) & (varDF.alt == i[1][2]) & varDF.feature == i[1][3]])
        frequency = float(float(count)/float(total))

def printRunInfo(config, outputDirName):
    '''
    Print useful information about the run, including the version, time, and configuration used
    '''
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
    TO DO: total runtime? or runtime breakdown per segment of the program?
    ofRunInfo.write("Variants Pre-filter: \n")
    ofRunInfo.write("        Post-filter: \n")
    '''
    ofRunInfo.close()
    return True

def printCounts2(outputDirName, varDF):
    '''
    Print counts per feature
    Replacement function
    '''

    # ============================================================
    # counts.txt
    try:
        ofCounts = open(outputDirName + "/counts.txt", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))

    grouped = varDF.groupby('feature')

    numHits = grouped['sample'].count()
    numHits.name = 'numHits'

    uniqueHits = grouped['sample'].count().groupby(level=0).count()
    uniqueHits.name = 'uniqueHits'

    numSamples = grouped['sample'].nunique()
    numSamples.name = 'numSamples'

    # TO DO: Construct DF from the 3 columns
    # TO DO: Write the DF to ofCounts
    
    return True

def printCounts(outputDirName, knownFeatures):
    '''
    Print counts per feature 
    '''
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

    ofCounts.close()
    print("\t{0}: {1} rows".format(ofCounts.name, nrow))
    return True

def collapseVariantDetails(group):
    '''
    Pandas operations to support the printVariantDetails family of functions. 
    Collapses variant rows that share the same contig, position, ref allele, alt allele, and feature.
    Input: a pandas groupby object
    Output: a pandas dataframe object
    '''
    outvals = []
    columns = list(group.keys().values)
    for column in columns:
        outstring = ''
        if column in ['vf', 'dp', 'sample', 'source']:  # the only columns that need to be concatenated 
                                                        # the others are uniquified and "always" yield 1 value
            for i in group[column].values:
                outstring += str(i) + ", "
            outvals.append( outstring[:-2] )            # trim the extra ', ' off the end of outstring 
        else:
            outvals.append( group[column].unique() )
    outDF = pd.DataFrame( dict(zip(columns,outvals)), columns=columns )
    return outDF

def printVariantDetails2(outputDirName, varDF, numSamples):
    '''
    Replacement function to print all information about each mutation,
    combining all mutations (irrespective of in how many samples they appear)
    into a single row
    Note: chrom, position, ref, alt, and feature are all required to uniquely identify a mutation 
          indels may have the same chr, pos, but different ref/alt
    '''

    # =========================================================
    # variant_details.txt
    try:
        ofVariantDetails = open(outputDirName + "/variant_details.txt", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))
    # Group by (chr, pos, ref, alt, feature)
    grouped = varDF.groupby(['chr', 'pos', 'ref', 'alt', 'feature'])
    # apply collapsing function to each pandas group
    out = grouped.apply(collapseVariantDetails)
    # print the new, collapsed dataframe to a file
    out.sort(['chr','pos','feature']).to_csv(ofVariantDetails, sep='\t', na_rep='?', index=False)

    return True

def printVariantDetails(outputDirName, knownFeatures, varDF, total):
    '''
    Print all information about each mutation. Mutations are combined into 1 row per location and ref/alt, with multiple samples in the 'Source' column when necessary.
    Output is in text format. 
    '''
    try:
        ofVariantDetails = open(outputDirName + "/variant_details.txt", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))

    # =========================================================
    # variant_details.txt

    ofVariantDetails.write('Feature\tContig\tPos\tRef\tAlt\tVF\tDP\t')
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
                    if len(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()) == 1:
                        ofVariantDetails.write(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()[0] + "\t") 
                    else:
                        tempdb = ""
                        for i in [x.strip(']').strip('[') for x in str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()).replace(',','').split(' ')]:
                            if str(i) not in str(tempdb):
                                tempdb += str(i) + ","
                        ofVariantDetails.write(tempdb.strip(',') + "\t")
                ofVariantDetails.write(str(len(var.source.split(','))) + "\t")
                ofVariantDetails.write(str(float(len(var.source.split(',')))/float(total)) + "\t" )
                ofVariantDetails.write(var.source + "\n")
                nrow += 1
                
    ofVariantDetails.close()
    print("\t{0}: {1} rows".format(ofVariantDetails.name, nrow))
    return True

def printLongVariantDetails(outputDirName, knownFeatures, varDF, total):
    '''
    Similar to printVariantDetails above, but writes each instance of a mutation to a new row. 
    Each mutation is written once per source instead of combining reoccurring mutations in to 1 unique row.
    '''
    try:
        ofLongVariantDetails = open(outputDirName + "/long_variant_details.txt", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))

    # =========================================================
    # long_variant_details.txt

    ofLongVariantDetails.write('Feature\tContig\tPos\tRef\tAlt\tVF\tDP\t')
    if SnpEff_switch:
        ofLongVariantDetails.write('Effect\tFC\t')
    if database_switch:
        ofLongVariantDetails.write('Annotation\t')
    ofLongVariantDetails.write('Count\tFrequency\tSource\n')

    masterList = list(knownFeatures.values())
    sortedList = sorted(masterList, key=lambda k: k.numVariants(), reverse=True)
    nrow = 0

    for feature in sortedList:
        if knownFeatures[feature.name].variants:
            for var in knownFeatures[feature.name].uniqueVariants():
                for j in range(len(var.dp.split(','))):
                    ofLongVariantDetails.write(feature.name + "\t")
                    ofLongVariantDetails.write(var.pos.chrom + "\t")
                    ofLongVariantDetails.write(str(var.pos.pos) + "\t")
                    ofLongVariantDetails.write(var.ref + "\t")
                    ofLongVariantDetails.write(var.alt + "\t")
                    ofLongVariantDetails.write(var.frac.split(', ')[j] + "\t")
                    ofLongVariantDetails.write(var.dp.split(', ')[j] + "\t")
                    
                    if SnpEff_switch:
                        if str(var.eff) != str(''):
                            ofLongVariantDetails.write(str([ x for x in set(str(var.eff).split(', '))]).replace("''","").replace(", ","").strip(']').strip('[').strip("'") + "\t")
                        else:
                            ofLongVariantDetails.write(str('?') + "\t")
                        if str(var.fc) != str(''):
                            ofLongVariantDetails.write(str([ x for x in set(str(var.fc).split(', '))]).replace("''","").replace(", ","").strip(']').strip('[').strip("'") + "\t")
                        else:
                            ofLongVariantDetails.write(str('?') + "\t")
                    if database_switch:
                        if len(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()) == 1:
                            ofLongVariantDetails.write(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()[0] + "\t") 
                        else:
                            tempdb = ""
                            for i in [x.strip(']').strip('[') for x in str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()).replace(',','').split(' ')]:
                                if str(i) not in str(tempdb):
                                    tempdb += str(i) + ","
                            ofLongVariantDetails.write(tempdb.strip(',') + "\t")
                    ofLongVariantDetails.write(str(len(var.source.split(','))) + "\t")
                    ofLongVariantDetails.write(str(float(len(var.source.split(',')))/float(total)) + "\t" )
                    ofLongVariantDetails.write(var.source.split(', ')[j] + "\n")
                    nrow += 1
                
    ofLongVariantDetails.close()
    print("\t{0}: {1} rows".format(ofLongVariantDetails.name, nrow))
    return True

def printVariantDetailsXLS(outputDirName, knownFeatures, varDF, total):
    '''
    Identical to the printVariantDetails function above, but in XLS format instead of TXT
    '''
    # =========================================================
    # variant_details.xls

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
                    if len(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()) == 1:
                        ofVariantDetails.write(nrow, ncol, str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()[0]) )
                        ncol += 1
                    else:
                        tempdb = ""
                        for i in [x.strip(']').strip('[') for x in str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()).replace(',','').split(' ')]:
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
                if nrow >= 65536:
                    throwWarning("Excel xls format cannot handle all your mutations")
                    throwWarning("{0} output could not be completed!".format(str(outputDirName + '/variant_details.xls')))
                    return
    workbook.save(outputDirName + '/variant_details.xls')
    print("\t{0}: {1} rows".format(str(outputDirName + '/variant_details.xls'), nrow))
    return True

def printVariantDetailsXLS2(outputDirName, varDF, numSamples):
    '''
    Replacement function to print all information about each mutation,
    combining all mutations (irrespective of in how many samples they appear)
    into a single row
    Note: chrom, position, ref, alt, and feature are all required to uniquely identify a mutation 
          indels may have the same chr, pos, but different ref/alt
    '''

    # =========================================================
    # variant_details.xls
    try:
        ofVariantDetails = ExcelWriter(str(outputDirName) + '/variant_details.xls')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))
    # Group by (chr, pos, ref, alt, feature)
    grouped = varDF.groupby(['chr', 'pos', 'ref', 'alt', 'feature'])
    # apply collapsing function to each pandas group
    out = grouped.apply(collapseVariantDetails)
    # print the new, collapsed dataframe to a file
    out.sort(['chr','pos','feature']).to_excel(ofVariantDetails, 'Variant Details', na_rep='?', index=False)
    ofVariantDetails.save()
    return True

def printLongVariantDetailsXLS2(varDF, outputDirName):
    ofLongVariantDetails = ExcelWriter(str(outputDirName) + '/long_variant_details.xls')
    varDF.sort(['chr','pos','feature']).to_excel(ofLongVariantDetails, 'Long Variant Details', na_rep='?', index=False)
    ofLongVariantDetails.save()

def printLongVariantDetailsXLS(outputDirName, knownFeatures, varDF, total):
    '''
    Similar to printVariantDetailsXLS above, but writes each instance of a mutation to a new row. 
    Each mutation is written once per source instead of combining reoccurring mutations in to 1 unique row.
    '''
    # =========================================================
    # long_variant_details.xls

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
                        if len(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()) == 1:
                            ofLongVariantDetails.write(nrow, ncol, str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()[0]) )
                            ncol += 1
                        else:
                            tempdb = ""
                            for i in [x.strip(']').strip('[') for x in str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()).replace(',','').split(' ')]:
                                if str(i) not in str(tempdb):
                                    tempdb += str(i) + ","
                            ofLongVariantDetails.write(nrow, ncol, str(tempdb.strip(',')))
                            ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, int(len(var.source.split(','))))
                    ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, float(float(len(set(var.source.split(', '))))/float(total))  )
                    ncol += 1
                    ofLongVariantDetails.write(nrow, ncol, str(var.source.split(', ')[j]))
                    ncol += 1
                    nrow += 1
                    if nrow >= 65536:
                        throwWarning("Excel xls format cannot handle all your mutations")
                        throwWarning("{0} output could not be completed!".format(str(outputDirName + '/long_variant_details.xls')))
                        return
                
    workbook.save(outputDirName + '/long_variant_details.xls')
    print("\t{0}: {1} rows".format(str(outputDirName + '/long_variant_details.xls'), nrow))
    return True

def printVariantBed(outputDirName, knownFeatures):
    '''
    Print bed file of the variant locations
    '''
    try:
        ofVariantBeds = open(outputDirName + "/variant_locations.bed", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))
    # =========================================================
    # variant_locations.bed
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
                ofVariantBeds.write("\n")
                nrow += 1
    ofVariantBeds.close()
    print("\t{0}: {1} rows".format(ofVariantBeds.name, nrow))
    return True

def getMetricsVCF(var, sourcefile):
    ''' 
    Return the depth and variant frequency for the given variant, according to the given source file 
    Does not necessarily have to have occurred in the source
        Returns 0,0 if it wasn't reported
    '''
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
    ''' 
    Print 1 vcf file as output
    There is 1 column per sample for the information they do not share (Ex. DP, VF)
    '''
    try:
        ofBigVCF = open(outputDirName + "/variant_calls.vcf", 'w+')
    except:
        abortWithMessage("Error opening output files in {0}/".format(outputDirName))
    # =========================================================
    # variant_calls.vcf
        
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
                    if len(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()) == 1:
                        ofVariantDetails.write(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()[0] + "\t") 
                    else:
                        tempdb = ""
                        for i in [x.strip(']').strip('[') for x in str(varDF[(varDF.pos == int(var.pos.pos)) & (varDF.chr == str(var.pos.chrom)) & (varDF.ref == str(var.ref))  & (varDF.alt == str(var.alt))].datab.unique()).replace(',','').split(' ')]:
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
    return True

def printOutput(config, outputDirName, knownFeatures, gas, varDF):
    '''Output statistics and variant details to the specified output directory.'''

    startTime = time.clock()
    print("\n=== Writing output files to {0}/ ===".format(outputDirName))

    total = len(set(config.samples))                # used in frequency calculation. config.samples will use sample count as the denominator. samples.inputFiles will use file count
    printRunInfo(config, outputDirName)
    
    if 'counts' in config.outputFormats:
        printCounts(outputDirName, knownFeatures)
    if 'txt' in config.outputFormats:
        printVariantDetails2(outputDirName, varDF, total)
    if 'bed' in config.outputFormats:
        print("bed format not yet supported")
        #printVariantBed(outputDirName, knownFeatures)
    if 'xlwt' in config.outputFormats and 'xlwt' in sys.modules:
        printVariantDetailsXLS(outputDirName, knownFeatures, varDF, total)
    if 'default' in config.outputFormats:
        printVariantDetails2(outputDirName, varDF, total)
        if 'xlwt' in sys.modules: printVariantDetailsXLS2(outputDirName, varDF, total)
        printCounts(outputDirName, knownFeatures)
    if 'long' in config.outputFormats and 'xlwt' in sys.modules:
        printLongVariantDetailsXLS2(varDF, outputDirName)
    if 'longtxt' in config.outputFormats:
        print("long text format not yet supported")
        #printLongVariantDetails(outputDirName, knownFeatures, varDF, total)
    if 'vcf' in config.outputFormats:
        print("long vcf format not yet supported")
        #printBigVCF(outputDirName, knownFeatures, varDF, config.inputFiles)

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

    varDF, knownFeatures, gas = parseVariantFiles(list(config.inputFiles), knownFeatures, gas, config.databases, config.filters, config.regions)

    printOutput(config, str(config.outputDir), knownFeatures, gas, varDF)
    
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
