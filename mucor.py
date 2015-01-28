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
import json

# nonstandard modules
import numpy as np
import pandas as pd
from pandas import ExcelWriter
import HTSeq

# optional modules
try:
    import tabix
except ImportError:
    print("Tabix module not found; database features disabled")

# mucor modules
import mucorfilters as mf
from variant import Variant
from mucorfeature import MucorFeature
import inputs
import output
from config import Config
from databases import dbLookup 
from info import Info

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
    '''
    Import the JSON config file from mucor_config.py. 
    Reads the config file into a dictionary, then writes each dictionary entry into the respective Config class position.
    '''
    config = Config()
    try:
        JD = json.load(open(json_config,'r'))
    except:
        abortWithMessage("Could not load the given JSON config file. See the example config for proper formatting.")

    # write dictionary values into a more human-friendly config class
    config.featureType = JD['feature']
    config.outputDir = JD['outputDir']
    config.union = JD['union']
    config.fast = JD['fast']
    config.gff = JD['gff']
    config.outputFormats = list(set(JD['outputFormats'])) # 'set' prevents repeated formats from being written multiple times
    if JD['databases']:
        config.databases = JD['databases']
    if str(JD['regions']):
        config.regions = JD['regions']
    else:
        config.regions = []
    # comma separated list of acceptable VCF filter column values
    config.filters = JD['filters']

    '''
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
    '''
    config.inputFiles = []
    config.samples = []
    for i in JD['samples']:
        config.samples.append(i['id'])
        for j in i['files']:
            filename = str(j['path']).split('/')[-1]
            config.filename2samples[filename] = i['id']
            config.source[filename] = j['source']
            '''
            if not SnpEff_switch and str(j['type']) == str('vcf') and bool(j['snpeff']) == bool(True):
                SnpEff_switch = bool(True)
            '''
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
    
    annotFileName = ".".join(gffFileName.split('/')[-1].split('.')[:-1])
    # pickled file is created for specific combinations of gff annotation and feature type. See below for more details. **
    archiveFilePath = str("/") + str(fast).strip('/') + str("/") + str(annotFileName) + str('_') + str(featureType) + str('.p')
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
    if bool(fast):
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
            if not os.path.exists( str("/".join(archiveFilePath.split('/')[:-1])) ):
                os.makedirs( str("/".join(archiveFilePath.split('/')[:-1])) )
            gas, knownFeatures, duplicateFeatures = constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures)
            archiveOut = open(archiveFilePath, 'wb')
            pickle.dump(gas, archiveOut, -1) ### pickle feature only works with full annotation files
            pickle.dump(knownFeatures, archiveOut, -1)
            pickle.dump(duplicateFeatures, archiveOut, -1)
            archiveOut.close()
    if not bool(fast):
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
        for kvar in knownFeatures[featureName].variants: # "kvar" stands for "known variant"
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

def groupCount(grp):
    ''' 
    function to count the number of samples that share the same mutation
    Input: Acts on a pandas groupby object
    Oputput: returns a DataFrame
    '''
    grp['count'] = len(grp['sample'])
    return grp

def parseVariantFiles(config, knownFeatures, gas, databases, filters, regions, total): 
    '''
    Read in all input files
    Record mutations in the variant dataframe 
    '''

    startTime = time.clock()

    # All variants stored in long (record) format
    # in a pandas dataframe
    #
    # However, it is MUCH faster to initially store them in a python Dict
    # Then convert to the pandas DF at the end
    varD = defaultdict(list) 
    variantFiles = list(config.inputFiles)

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

        while str(row).split("'")[1][0:2] == '##':
            row = varReader.next()

        header = row
        if len(header) == 0: raise ValueError('Invalid header')
        fieldId = dict(zip(header, range(0, len(header))))

        # after reading the two header rows, read data
        for row in itertools.islice(varReader, None):
            if filterRow(row, fieldId, filters, kind):  # filter rows as they come in, to prevent them from entering the dataframe
                continue                                # this allows us to print the dataframe directly and have consistent output with variant_details.txt, etc.

            # attempt to extract 'effect' and 'functional consequence' from the VCF line
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
                # reformat the lists to exclude blanks and be semicolon delimited
                for mut in set(muts):
                    if str(mut) != "":
                        effect += str(mut) + ";"
                for loc in set(loca):
                    if str(loc) != "":
                        fc += str(loc) + ";"
            except KeyError:
                # this VCF may not be snpEff annotated
                pass

            parser = inputs.Parser()
            source = config.source[fn.split('/')[-1]]
            var = parser.parse(source, row, fieldId, header, fn, effect, fc)
            if regions and not inRegionDict(var.pos.chrom, int(var.pos.pos), int(var.pos.pos), regionDict ):
                continue

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
            sample = config.filename2samples[str(fn.split('/')[-1])]
            source = os.path.split(fn)[1]
            dbEntries = dbLookup(var, databases)
            count = 0 # initialize this database column now to save for later
            freq = 0.0 # initialize this database column now to save for later

            # build dict to insert
            # Step 1: define the columns and their values
            # Step 1b.The database columns are variable, ranging from zero to many
            #         and must be inserted dynamically
            # Step 2. zip the column names and column values together into a dictionary
            # Step 3. Add this round to the master variants data frame dictinoary

            for feature in features.split(', '):
                # Removing columns from the following 'columns' list will mask them from output
                columns = ['chr','pos','ref','alt','vf','dp','feature','effect','fc']
                values  = [ chr,  pos,  ref,  alt,  vf,  dp,  feature,  effect,  fc]

                for dbName in sorted(dbEntries.keys()):
                    columns.append(dbName)
                    values.append(dbEntries[dbName])
                columns += ['count', 'freq', 'sample','source']
                values  += [ count,   freq,   sample,  source ]

                vardata = dict(zip( columns, values ))
                for key in vardata.keys():
                    varD[key].append(vardata[key])

        totalTime = time.clock() - startTime
        print("{0:02d}:{1:02d}\t{2}".format(int(totalTime/60), int(totalTime % 60), fn))
    
    # Transform data frame dictionary into pandas DF. Major speed increase relative to appending the DF once per variant
    varDF = pd.DataFrame(varD, columns=columns)
    # Dataframe operation to count the number of samples that exhibit each mutation
    varDF = varDF.groupby(['chr', 'pos', 'ref', 'alt', 'feature']).apply(groupCount)
    # Divide the count for each row in varDF by the total sample count from the JSON config file. Represents the percent of the input samples that exhibit each mutation.
    freqdict = {}
    for i in set(varDF['count'].values):
        freqdict[i] = float(i)/float(total)
    varDF['freq'] = varDF['count'].map(freqdict)
    # Clean up variant dataframe a little
    # position should be integer, not float
    varDF.pos   = varDF.pos.astype(int)
    # blank features should be more descriptive. this replaces empty strings with 'NO_FEATURE'. 
    #   this is important to later functions that put feature in a pandas index. na_rep pandas function only changes dataframe values, not the index.
    #   see output.FeatureXSample() for an example of when 'feature' is used in a pandas index.
    varDF['feature'][varDF.feature == ''] = 'NO_FEATURE'

    return varDF, knownFeatures, gas 

def printOutput(config, outputDirName, varDF):
    '''Output statistics and variant details to the specified output directory.'''

    startTime = time.clock()
    print("\n=== Writing output files to {0}/ ===".format(outputDirName))
    ow = output.Writer()
    for format in config.outputFormats:
        ow.write(varDF,format,outputDirName,config)

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

    # Total is used in frequency calculations. 
    #   supplying config.samples will use sample count as the denominator [canonical operation, ie: comparing samples]
    #   or, using samples.inputFiles will use file count [non-canonical operation, ie: comparing tools, or otherwise having many vcf files and 1 sample ID]
    total = len(set(config.samples))

    knownFeatures, gas = parseGffFile(str(config.gff), str(config.featureType), config.fast)

    varDF, knownFeatures, gas = parseVariantFiles(config, knownFeatures, gas, config.databases, config.filters, config.regions, total)

    printOutput(config, str(config.outputDir), varDF)
    
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
    
    # pretty print newline before exit
    print()


if __name__ == "__main__":
    if sys.hexversion < 0x02070000:
        raise RuntimeWarning("mucor should be run on python 2.7.0 or greater.")
    main()
