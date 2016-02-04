#!/usr/bin/env python
# -*- coding: utf8
#    Copyright 2013-2015 James S Blachly, MD and The Ohio State University
#
#    This file is part of Mucor.
#
#    Mucor is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Mucor is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Mucor.  If not, see <http://www.gnu.org/licenses/>.

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
import json

# nonstandard, required modules
try:
    import numpy as np
except ImportError as err:
    print(str(err) + ". This module is required.")
    sys.exit()
try:
    import pandas as pd
except ImportError as err:
    print(str(err) + ". This module is required.")
    sys.exit()
try:
    import HTSeq
except ImportError as err:
    print(str(err) + ". This module is required.")
    print("* NOTE! If you installed with pip install mucor,")
    print("* the numpy and pandas dependencies were automatically")
    print("* installed. However, due to a bug in the HTSeq")
    print("* installation, it cannot be automatically installed")
    print("* at the same time as numpy.")
    print("* SOLUTION: pip install HTSeq")
    sys.exit()

# optional modules
# will throw warnings later in the program 
try:
    import tabix
except ImportError:
    pass

try:
    import xlsxwriter
except ImportError:
    pass

# mucor modules
import mucorfilters as mf
from variant import Variant
from mucorfeature import MucorFeature
import inputs
import output
from config import Config
from databases import dbLookup,checkAndOpen 
from info import Info

def abortWithMessage(message):
    print("*** FATAL ERROR: " + message + " ***")
    exit(2)

def throwWarning(message, help = False):
    print("*** WARNING: " + message + " ***")
    return

def constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures, union):
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

        # skip full transcripts, full genes, and other items that are not exon ranges
        # including these lines will create '--union'-like behavior, regardless of whether the user passed the union option
        if feature.type != 'exon':
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
                if union:
                    # replace the start and end coordinates when adding SUCCESSIVE bits of a feature (e.g. exons)
                    # results in one, large feature, starting with the earliest upsream location and ending with the latest downstream location
                    if knownFeatures[feat.name].iv.start < feat.iv.start:
                        feat.iv.start = knownFeatures[feat.name].iv.start
                    if knownFeatures[feat.name].iv.end > feat.iv.end:
                        feat.iv.end = knownFeatures[feat.name].iv.end
                else:
                    pass # no-union - this does overwrite previous coordinates in knownFeatures,
                         # but should not matter as the actual coordinates are obtaind from 'gas'.
                         # Locations held in knownFeatures should not be used to identify the start and stop of a feature,
                         # since these locations will only reveal the last known region for the feature. 
                         # Later in the program, querying knownFeatures for the variants in a feature (ie: skipThisIndel function)
                         # may return variants that appear to be located in positions outside of the feature region reported by knownFeatures. 
                         # This is also caused by the 'no-union' overwrite, since the knownFeature region locations do not represent the whole feature.

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
    except ValueError as json_error:
        throwWarning(json_error.message)
        abortWithMessage("Could not load the given JSON config file. See the example config for proper formatting.")

    # write dictionary values into a more human-friendly config class
    config.featureType = JD['feature']
    config.outputDir = os.path.expanduser(JD['outputDir'])
    config.union = JD['union']
    config.fast = JD['fast']
    config.gff = JD['gff']
    config.outputFormats = list(set(JD['outputFormats'])) # 'set' prevents repeated formats from being written multiple times
    if JD['databases']:
        if 'tabix' in sys.modules: # make sure tabix is imported 
            for name,db in JD['databases'].items():
                dbPointer = checkAndOpen(db)
                if dbPointer:
                    #check for non-null pointers
                    config.databases[name] = dbPointer
        else: # the user supplied databases but did not successfully import tabix
            throwWarning("tabix module not found; database features disabled")
    if str(JD['regions']):
        config.regions = JD['regions']
    else:
        config.regions = []
    # comma separated list of acceptable VCF filter column values
    config.filters = JD['filters']

    config.inputFiles = []
    config.samples = []
    for i in JD['samples']:
        config.samples.append(i['id'])
        for j in i['files']:
            if j['type'] == "bam":
                # skip the bam files
                continue
            filename = os.path.basename(j['path'])
            config.filename2samples[filename] = i['id']
            config.source[filename] = j['source']
            config.inputFiles.append(j['path'])

    return config 

def parseGffFile(gffFileName, featureType, fast, union):
    '''
    Parse the GFF/GTF file. Return tuple (knownFeatures, GenomicArrayOfSets)
    Haplotype contigs are explicitly excluded because of a coordinate crash (begin > end)
    '''
    
    # TO DO: command line flag should indicate that variants in INTRONS are counted
    # This is called --union, see below
    
    startTime = time.clock()
    print("\n=== Reading GFF/GTF file {0} ===".format(gffFileName))
    
    if union:
        unionstatus = str("union")
    else:
        unionstatus = str("no_union")
    annotFileName = os.path.splitext(os.path.basename(gffFileName))[0]
    # pickled file is created for specific combinations of gff annotation and feature type. See below for more details. **

    archiveFilePath = os.path.expanduser(fast) + str("/") + str(annotFileName) + str('_') + str(featureType) + str("_") + unionstatus + str('.p')
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
    # ** These three items will change depending on the supplied gff annotation AND the feature selected AND whether union was used. 
    #    Thus, each pickle file will be for a specific combination of gff, feature, and union status, hence the naming convention.
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
            # no pickled annotation exists for this combination of gff and feature; creating it in the directory provided to the 'fast' option
            print("Cannot locate annotation archive for " + os.path.basename(gffFileName) + str(" w/ ") + str(featureType) + " and --union=" + str(union) )
            print("   Reading in annotation and saving archive for faster future runs") 
            if not os.path.exists( str("/".join(archiveFilePath.split('/')[:-1])) ):
                os.makedirs( str("/".join(archiveFilePath.split('/')[:-1])) )
            gas, knownFeatures, duplicateFeatures = constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures, union)
            archiveOut = open(archiveFilePath, 'wb')
            pickle.dump(gas, archiveOut, -1)
            pickle.dump(knownFeatures, archiveOut, -1)
            pickle.dump(duplicateFeatures, archiveOut, -1)
            archiveOut.close()
    if not bool(fast):
    # ignore pickles function entirely. Won't check for it and won't attempt to create it
        gas, knownFeatures, duplicateFeatures = constructGAS(gffFile, featureType, knownFeatures, duplicateFeatures, union)

    if duplicateFeatures:
        print("*** WARNING: {0} {1}s found on more than one contig".format(len(duplicateFeatures), featureType))
    # remove gene-centric feature names that are known to have problematic bins caused by multiple copies on the same gene
    # run the accompanying 'detect union bin errors' python script to generate a list of incompatible genes
    if union and "gene" in featureType:
        try:
            badGenes = open('union_incompatible_genes.txt')
            badGeneList = set()
            print("Removing genes with multiple copies on the same contig, which cause incorrect feature bins with 'union'\n\tUsing 'union_incompatible_genes.txt'")
            for line in badGenes:
                badGene = line.strip()
                badGeneList.add(badGene)
            for interval, sets in gas.steps():
                gas[interval] = set( [ x for x in sets if x not in badGeneList ] )
            for badGene in badGeneList:
                try:
                    del knownFeatures[badGene]
                    # uncomment below to see which genes are being deleted 
                    # print(badGene + " deleted from knownfeat")
                except KeyError:
                    pass
        except IOError:
            # There likely isn't a list of genes
            # Recommend that the user obtain said list, or face potential errors in feature labels and counts
            pass
    else:
        # not using a gene_name-type of feature in combination with --union
        # no known issues should arise with this configuration 
        pass

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
    if kind in ["vcf", "vcf.gz"]:
        # this is handled in 1 line elsewhere in the program. Specifically, within the parseVariantFiles function, under the [vcf, vcf.gz] block
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
    
    for kvar in knownFeatures[featureName].variants: # "kvar" stands for "known variant"
        if kvar.pos.pos == var.pos.pos and kvar.pos.chrom == var.pos.chrom and \
        ( kvar.ref != var.ref or kvar.alt != var.alt ) and \
        indelDelta(var.ref,var.alt)[0] == indelDelta(kvar.ref, kvar.alt)[0] and \
        indelDelta(var.ref,var.alt)[1] == indelDelta(kvar.ref, kvar.alt)[1]:
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
    DEPRICATED: This has been reduced to a faster lambda function 
    '''
    return len(grp['sample'].unique())

def annotateDF(grp, databases):
    '''
    function to annotate variant dataframe with user-supplied vcf databases
    '''
    chrom = grp['chr'].unique()[0]
    pos = grp['pos'].unique()[0]
    ref = grp['ref'].unique()[0]
    alt = grp['alt'].unique()[0]
    var = Variant(source=None, sample=None, pos=HTSeq.GenomicPosition(chrom, pos), ref=ref, alt=alt, frac=None, dp=None, eff=None, fc=None)
    dbEntries = dbLookup(var, databases)
    return pd.Series(dbEntries)

def integrateVar(var, varD, config, gas, knownFeatures, unrecognizedContigs, unrecognizedMutations):
    '''
    Ingest the given variant object into the variant dictionary 
    '''
    # find bin for variant location
    try:
        resultSet = gas[ var.pos ]      # returns a set of zero to n IDs (e.g. gene symbols)
    except KeyError:                    # which I'll use as a key on the knownFeatures dict
                                        # and each feature with matching ID gets allocated the variant
        # this mutation is on a contig unknown to the GAS
        resultSet = set()
        unrecognizedContigs.add(var.pos.chrom)
        unrecognizedMutations += 1
   
    if resultSet:
        for featureName in resultSet:
            if len(var.ref) != len(var.alt): # confirm that this mutation is an indel, not a snp
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
    features = ', '.join( resultSet )   # join with comma to handle overlapping features
    effect = var.eff
    fc = var.fc
    sample = var.sample
    source = var.source
    count = 0 # initialize this database column now to save for later
    freq = 0.0 # initialize this database column now to save for later
    
    # build dict to insert
    # Step 1: define the columns and their values
    # Step 2. zip the column names and column values together into a dictionary
    # Step 3. Add this round to the master variants data frame dictinoary

    for feature in features.split(', '):
        # Removing columns from the following 'columns' list will mask them from output
        columns = ['chr','pos','ref','alt','vf','dp','feature','effect','fc','count','freq','sample','source']
        values  = [ chr,  pos,  ref,  alt,  vf,  dp,  feature,  effect,  fc , count,  freq,  sample, source  ]

        vardata = dict(zip( columns, values ))
        for key in vardata.keys():
            varD[key].append(vardata[key])
    return varD, unrecognizedContigs, unrecognizedMutations

def parseVariantFiles(config, knownFeatures, gas, databases, filters, regions, total) :
    '''
    Read in all input files
    Record mutations in the variant dataframe 
    '''

    startTime = time.clock()
    unrecognizedContigs = set()
    unrecognizedMutations = 0

    # All variants stored in long (record) format
    # in a pandas dataframe
    #
    # However, it is MUCH faster to initially store them in a python Dict
    # Then convert to the pandas DF at the end
    varD = defaultdict(list) 
    variantFiles = set(list(config.inputFiles)) # uniquify the file list, in the case of the same multi-sample VCF being defined for multiple samples

    if regions: # has the user specified any particular regions or region files to focus on?
        regionDict = defaultdict(set)
        for item in regions:
            if os.path.splitext(item)[1].lower() == 'bed':      # this item is a bed file
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
        kind = str( os.path.splitext(fn)[-1].strip('.').lower() )
        if kind == "gz" and fn.endswith(".vcf.gz"):
            kind = "vcf.gz"
        try:
            varFile = open(fn, 'rb')
        except IOError:
            # file does not exist or cannot be opened
            throwWarning('{} could not be opened'.format(fn))
            continue    # next fn in variantFiles

        if kind == "out":
            # parse as mutect '.out' type format
            varReader = csv.reader(varFile, delimiter="\t")
            row = varReader.next()
            while row[0].startswith('##'):
                row = varReader.next()
            header = row
            if len(header) == 0: raise ValueError('Invalid header')
            fieldId = dict(zip(header, range(0, len(header))))
            for row in itertools.islice(varReader, None):
                if filterRow(row, fieldId, filters, kind):  # filter rows as they come in, to prevent them from entering the dataframe
                    continue                                # this allows us to print the dataframe directly and have consistent output with variant_details.txt, etc.
                source = config.source[ os.path.basename(fn) ]
                parser = inputs.Parser()
                parser.row = row
                parser.source = source
                parser.fieldId = fieldId
                parser.header = header
                parser.fn = fn 
                var = parser.parse_MuTectOUT()
                var.sample = config.filename2samples[os.path.basename(fn)]
                if not var:
                    # this mutation had no data in this sample
                    continue
                if regions and not inRegionDict(var.pos.chrom, int(var.pos.pos), int(var.pos.pos), regionDict ):
                    continue
                varD, unrecognizedContigs, unrecognizedMutations = integrateVar(var, varD, config, gas, knownFeatures, unrecognizedContigs, unrecognizedMutations)

        elif kind in ["vcf", "vcf.gz"]:
            # start htseq vcf reader
            varReader = HTSeq.VCF_Reader(str(fn))
            varReader.parse_meta()
            varReader.make_info_dict()
            for row in varReader:
                if row.filter not in filters:
                    continue
                if regions and not inRegionDict(row.pos.chrom, int(row.pos.pos), int(row.pos.pos), regionDict ):
                    continue
                row.unpack_info(varReader.infodict)
                parser = inputs.Parser()
                source = config.source[ os.path.basename(fn) ]
                try:
                    samps = parser.parse(row, source)
                except KeyError:
                    throwWarning("parsing " + fn + " as '" + source + "'")
                    print("Cannot parse file from an unsupported or unknown variant caller. \nPlease use supported variant software, or compose an input module compatible with inputs.py")
                    print("Options include: " + str(parser.supported_formats.keys()).replace("'",""))
                    break
                for var in samps:
                    if len(samps) == 1:
                        # if not multi-sample VCF, pull sample ID from the json config
                        var.sample = config.filename2samples[ os.path.basename(fn) ]
                    elif len(samps) > 1:
                        # if multi-sample VCF, use sample ID as defined by the VCF column(s) rather than the config
                        pass
                    if not var or var.sample not in config.samples:
                        # this sample has no mutation data at the given location, or this sample was not specified as a sample of interest in the JSON config 
                        continue
                    var.source = os.path.basename(fn)
                    
                    varD, unrecognizedContigs, unrecognizedMutations = integrateVar(var, varD, config, gas, knownFeatures, unrecognizedContigs, unrecognizedMutations)
        else:
            throwWarning("Unable to parse file with extension '{0}': {1}".format(kind, fn))
            continue
        if unrecognizedContigs:
            throwWarning("{0} Contigs and {1} mutations are in areas unknown to the genomic array of sets. If using --archive_directory, perhaps try again without it.".format( len(unrecognizedContigs), unrecognizedMutations ))
        totalTime = time.clock() - startTime
        print("{0:02d}:{1:02d}\t{2}".format(int(totalTime/60), int(totalTime % 60), fn))
    
    columns = ['chr','pos','ref','alt','vf','dp','feature','effect','fc','count','freq','sample','source']
    # Transform data frame dictionary into pandas DF. Major speed increase compared to appending variants to the DF while reading the input files. 
    try:
        varDF = pd.DataFrame(varD, columns=columns)
        # Drop samples with no detected mutation at the given loc
        varDF = varDF.dropna(subset=['dp','vf'], how='all')
        # Cast dp, vf, and pos as appropriate types
        varDF['dp'] = varDF['dp'].astype(int)
        varDF['vf'] = varDF['vf'].astype(float)
        varDF.pos   = varDF.pos.astype(int)
    except UnboundLocalError:
        if not varD:            
            # variant dictionary is empty, implying no variants were encountered and kept
            abortWithMessage("Variant DataFrame is empty! No mutations passed filter; check FILTER column of input VCF against allowed filters in JSON config.")
        else:
            abortWithMessage("Could not convert Variant Dictionary into a DataFrame")
    except ValueError:
        # likely couldn't convert NaN values from dp column to int
        # proceed anyway, treating them as floats instead
        varDF['dp'] = varDF['dp'].astype(float)
    if len(varDF) == 0:
        abortWithMessage("Variant DataFrame is empty! No mutations passed filter; check FILTER column of input VCF against allowed filters in JSON config.")
    # Annotate dataframe using user-supplied vcf databases
    if databases:
        print("\n=== Comparing Your Variants to Known VCF Databases ===")
        startTime = time.clock()
        annotSeries = varDF.groupby(['chr', 'pos', 'ref', 'alt']).apply( annotateDF, databases=databases )
        annotDF = pd.DataFrame(annotSeries)
        new_index = pd.Index([x for x in varDF.columns[:-4]] + [x for x in annotDF.columns] + [x for x in varDF.columns[-4:]])
        varDF = varDF.merge(annotDF, left_on=["chr", "pos", "ref", "alt"], right_index=True)[new_index]
        totalTime = time.clock() - startTime
        print("{0:02d}:{1:02d}".format(int(totalTime/60), int(totalTime % 60)))

    # Delete VAF columns from databases that do not have population variant allele frequencies recorded
    # These will manifest as a whole VAF column full of '?' 
    for column in varDF.keys():
        if "_VAF" in str(column) and len(set(varDF[column])) == 1:
            varDF.drop(column, axis=1, inplace=True)
 
    # Dataframe operation to count the number of samples that exhibit each mutation
    countSeries = varDF.groupby(['chr', 'pos', 'alt']).apply(lambda x: len(set(x['sample'])))
    countDF = pd.DataFrame(countSeries, columns=['count'])
    # This appears messy, but runs very fast regardless. It maps the grouped count dataframe back to the variant dataframe 
    varDF['count'] = varDF[['chr','pos','alt']].merge(countDF, left_on=['chr','pos','alt'], right_index=True)['count']
    # Divide the count for each row in varDF by the total sample count from the JSON config file. Represents the percent of the input samples that exhibit each mutation.
    freqdict = {}
    for i in set(varDF['count'].values):
        freqdict[i] = float(i)/float(total)
    varDF['freq'] = varDF['count'].map(freqdict)

    # blank features should be more descriptive. this replaces empty strings with 'NO_FEATURE'. 
    #   this is important to later functions that put feature in a pandas index. na_rep pandas function only changes dataframe values, not the index.
    #   see output.FeatureXSample() for an example of when 'feature' is used in a pandas index.
    
    # Avoid chained indexing:
    # use DataFrame.loc[row_index, col_index]
    varDF.loc[ varDF.feature == '', 'feature'] = 'NO_FEATURE'
    
    # replace all remaining blank cells with '?'

    #stop() # this command throws a warning
    varDF.replace('', '?', inplace=True)

    return varDF, knownFeatures, gas 

def printOutput(config, outputDirName, varDF):
    '''Output run statistics and variant details to the specified output directory.'''

    startTime = time.clock()
    print("\n=== Writing output files to {0}/ ===".format(outputDirName))
    ow = output.Writer()
    # identify formats ending in 'xlsx' as the "excel formats," requiring XlsxWriter
    excel_formats = [ plugin_name for plugin_name,file_name in ow.file_names.items() if os.path.splitext(file_name)[1]==".xlsx" ]
    if "all" in config.outputFormats:
        # replace output type 'all' with a lsit of all supported output types
        #    and remove 'all' and 'default' to prevent recursive execution of the modules
        config.outputFormats = ow.supported_formats.keys()
        config.outputFormats.remove('all')
        config.outputFormats.remove('default')
        if 'xlsxwriter' not in sys.modules and excel_formats:
            # if xlsxwriter is not present and the user selected excel output formats, remove excel formats from output formats
            config.outputFormats = [ x for x in config.outputFormats if x not in excel_formats ]
            throwWarning("xlsxwriter module not found; Excel outputs disabled")

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

    parser = argparse.ArgumentParser()
    parser.add_argument("json", help="Pass 1 json config as an argument")
    args = parser.parse_args()
    config = parseJSON(args.json)

    if not os.path.exists(config.gff):
        abortWithMessage("Could not find GFF file {0}".format(config.gff))
    if os.path.exists(config.outputDir) and [x for x in os.listdir(config.outputDir) if x in output.Writer().file_names.values() ]:
        abortWithMessage("The directory {0} already exists and contains output. Will not overwrite.".format(config.outputDir))
    elif not os.path.exists(config.outputDir):
        try:
            os.makedirs(config.outputDir)
        except OSError:
            abortWithMessage("Error when creating output directory {0}".format(config.outputDir))

    # check that all specified variant files exist
    for fn in config.inputFiles:
        if not os.path.exists(fn):
            abortWithMessage("Could not find variant file: {0}".format(fn))

    # Total is used in frequency calculations. 
    #   supplying config.samples will use sample count as the denominator [canonical operation, ie: comparing samples]
    #   or, using samples.inputFiles will use file count [non-canonical operation, ie: comparing tools, or otherwise comparing many vcf files with no regard for sample ID]
    total = len(set(config.samples))

    knownFeatures, gas = parseGffFile(str(config.gff), str(config.featureType), config.fast, config.union)
    varDF, knownFeatures, gas = parseVariantFiles(config, knownFeatures, gas, config.databases, config.filters, config.regions, total)
    printOutput(config, str(config.outputDir), varDF)
    
    # pretty print newline before exit
    print()


if __name__ == "__main__":
    if sys.hexversion < 0x02070000:
        raise RuntimeWarning("mucor should be run on python 2.7.0 or greater.")
    main()
