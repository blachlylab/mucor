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

# let print() be a function rather than statement
# ala python3
from __future__ import print_function

# python standard modules
import os
import re
import sys
reload(sys)  
sys.setdefaultencoding('utf8')
import time
import argparse
import csv
import itertools
from collections import defaultdict
import gzip
import cPickle as pickle
import json

# nonstandard modules
import numpy as np
import pandas as pd
from pandas import ExcelWriter
import HTSeq
import pysam

# optional modules

# mucor modules
from config import Config
from mucor_config import parseAndValidateRegions 

def abortWithMessage(message):
    print("*** FATAL ERROR: " + message + " ***")
    exit(2)

def throwWarning(message, help = False):
    print("*** WARNING: " + message + " ***")
    return

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
    config.outputDir = os.path.expanduser(JD['outputDir'])
    
    # no longer take region definition from json config file
    #   command line argument now
    if str(JD['regions']):
        config.regions = JD['regions']
    else:
        config.regions = []

    config.inputFiles = []
    config.samples = []
    config.bams = {}
    for i in JD['samples']:
        # make a list of bam file paths defined by the JSON config
        #   i.e. extract the path of every file, if the file type is bam
        bams = [x['path'] for x in i['files'] if x['type'] == "bam"]
        if len(bams) >= 1:
            for bam in bams:
                config.inputFiles.append(bam)
                config.samples.append(i['id'])
                try:
                    config.bams[i['id']].append(bam)
                except KeyError:
                    config.bams[i['id']] = []
                    config.bams[i['id']].append(bam)
        else:
            throwWarning("Sample {0} has no BAM file!".format(i['id']))
    if not config.bams:
        abortWithMessage("No bam files defined!")
    return config 

def parseRegionBed(regionfile, regionDictIn):
    ''' 
    Read through the supplied bed file and add the rows to the input region dictionary before returning it. Appends to the input dict, rather than overwriting it.
    '''
    regionDict = regionDictIn
    for line in open(str(regionfile),'r'):
        col = line.strip().split("\t")
        chrom = col[0]
        start = col[1]
        end = col[2]
        if len(col) > 3:
            name = col[3]
        else:
            name = str(chrom) + ":" + str(start) + "-" + str(end)
        regionDict[chrom].add((start,end,name))
    return regionDict

def inRegionDict(chrom, start, end, regionDict):
    '''Checks the given mutation location to see if it is in the dictionary of regions'''
    if regionDict[chrom]: # are there any regions of interest in the same chromosome as this mutation?
        for locs in regionDict[chrom]: # breaking the list of regions according to chromosome should significantly decrease the number of comparisons necessary 
            if locs[0] == 0 and locs[1] == 0: # chrN:0-0 is used to define an entire chromosome as a region of interest. 
                return True
            elif int(start) >= int(locs[0]) and int(end) <= int(locs[1]):
                return True
    return False

def GaugeDepth(config) :
    '''
    assess the depth of each sample at the given regions
    '''

    startTime = time.clock()
    regionDict = defaultdict(set)
    for item in config.regions:
        if str(str(item).split('.')[-1]).lower() == 'bed':      # this item is a bed file
            regionDict = parseRegionBed(item, regionDict)
        else:                                                   # this was defined as a string
            reg_chr = item[0]
            reg_str = item[1]
            reg_end  = item[2]
            if reg_str >= 0 and reg_end > reg_str:
                name = str(reg_chr) + ":" + str(reg_str) + "-" + str(reg_end)
            elif reg_str >= 0 and reg_str == reg_end:
                name = str(reg_chr) + ":" + str(reg_str)
            elif not reg_str and not reg_end:
                name = str(reg_chr)
            regionDict[reg_chr].add((reg_str, reg_end, name))

    if not regionDict:
        abortWithMessage("Regions not set!")

    covD = {'chr':[], 'start':[], 'stop':[], 'name':[], 'sample':[],'depth':[]}
    print("\n=== Reading BAM Files ===")
    for sid, fns in config.bams.items():
        # loop over all samples
        for fn in fns:
            # loop over all bam files for this sample
            try:
                samfile = pysam.AlignmentFile(fn, "rb" )
            except ValueError:
                throwWarning("Cannot open file {0}".format(fn))
                continue
            contig_length_dict = dict(zip(samfile.references, samfile.lengths)) # save contig lengths for later
            for contig, ROI in regionDict.items():
                for window in ROI:
                    bed_name = window[2]
                    # make window 0-based
                    try:
                        window = [int(window[0]) - 1, int(window[1])]
                    except TypeError: 
                        # start and/or stop were not defined - pull them from contig dictionary
                        if not window[0]:
                            start = 0
                        else:
                            start = int(window[0] - 1)
                        if not window[1]:
                            end = contig_length_dict[contig]
                        else:
                            end = int(window[1])
                        window = [start, end]
                    # loop over all ROIs, checking this bam 
                    if config.p:
                        #point method 
                        tmp_dict = {}
                        position = round((window[1] - window[0])/2.0) + window[0]
                        avg_covg = samfile.count(contig, position - 1, position)

                        #for position in range(window[0],window[1]):
                        #    region = str(contig) + ':' + str(position) + '-' + str(position)
                        #    tmp_dict[position] = samfile.count(region=region)
                        #avg_covg = np.mean(tmp_dict.values())
                    elif config.c:
                        #read count method
                        avg_covg = samfile.count(contig, window[0], window[1])
                        #note that "avg_covg" is only a name here - it is the total count of reads, not an average! 
                    else:
                        #complete average method
                        #'''
                        tmp_dict = {}
                        for position in range(window[0],window[1]):
                            tmp_dict[position] = 0
                        for pileupcolumn in samfile.pileup(contig, window[0], window[1],stepper='all'):
                            #loop over reads that hit the window locations and record coverage 
                            # 'stepper = all' yields mapped, primary, non-duplicate (identified by sam flag), QC pass reads
                            try:
                                tmp_dict[pileupcolumn.pos]
                                tmp_dict[pileupcolumn.pos] = pileupcolumn.n       
                            except KeyError:
                                #skip this position if it's not in the region dict
                                continue
                        avg_covg = np.mean(tmp_dict.values())
                        #'''
                        '''
                        #this behaves erratically and does not produce the same number if run repeatedly
                        #   not sure how to use the function, but it could be faster than pileup
                        counter = 0
                        for ct_cov in samfile.count_coverage(contig, window[0], window[1], read_callback = 'all'):
                            for nt_arr in ct_cov:
                                counter += int(nt_arr)
                        stop()
                        avg_covg = counter/float(window[1] - window[0])
                        '''


                    covD['chr'].append(str(contig))
                    covD['start'].append(int(window[0]) + 1)
                    covD['stop'].append(int(window[1]))
                    covD['name'].append(str(bed_name))
                    covD['sample'].append(str(sid))
                    covD['depth'].append(float(avg_covg))     
            samfile.close()
            totalTime = time.clock() - startTime
            print("{0:02d}:{1:02d}\t{2}".format(int(totalTime/60), int(totalTime % 60), fn))
    covDF = pd.DataFrame.from_dict(covD)[['chr','start','stop','name','sample','depth']]
    covDF = covDF.groupby(['chr','start','stop','name','sample'])['depth'].apply(sum).reset_index()
            
    totalTime = time.clock() - startTime
    print("\n{0:02d}:{1:02d}\t{2}".format(int(totalTime/60), int(totalTime % 60), "Done"))
    return covDF

def SamplesToColumns(grp):
    '''helper function to rearrange the coverage dataframe'''
    columns = ['chr','start','stop','name']
    values = [grp['chr'].values[0], grp['start'].values[0], grp['stop'].values[0],grp['name'].values[0]]

    grp_dict = {}
    for i in grp[['sample','depth']].T.itertuples(): 
        grp_dict[i[0]] = list(i[1:])
    columns += grp_dict['sample']
    values += grp_dict['depth']
    return pd.DataFrame.from_dict(dict(zip(columns,values), index=[0]))[columns]


def printOutput(config, covDF):
    '''reformat the coverage dataframe and write it to xlsx file'''

    startTime = time.clock()
    print("\n=== Writing output files to {0}/ ===".format(config.outputDir))
    outDF = covDF.groupby(['chr','start','stop','name']).apply(SamplesToColumns).reset_index(drop=True)
    # pandas can actually start up an excel writer object using a non-existant or unwritable file path. The error will come when saving to such a path. tested in pandas version 0.16.2
    outXLSX = pd.ExcelWriter(str(config.outputDir) + "/" + config.outfile, engine='xlsxwriter')
    outDF.to_excel(outXLSX, 'Depth Gauge', index=False)
    try:
        outXLSX.save()
    except IOError:
        abortWithMessage("Output directory is not writable or doesn't exist: {0}".format(config.outputDir))
    totalTime = time.clock() - startTime
    print("\t{0}: {1} rows".format(str(config.outputDir) + "/" + config.outfile , len(outDF)))     
    return True

def main():
    print("\n=== Run Info ===")
    print("\t{0}".format(time.ctime() ) )
    print()

    parser = argparse.ArgumentParser()
    parser.add_argument("json", help="Pass 1 json config as an argument")
    parser.add_argument("-p", "--point", dest='p', action="store_true", help="point-location coverage approximation on middle of region windows")
    parser.add_argument("-c", "--count", dest='c', action="store_true", help="report total number of reads in each region window (not an average)")
    parser.add_argument("-r", "--regions", default="", help="Comma separated list of bed regions and/or bed files by which to limit output. Ex: chr1:10230-10240,chr2,my_regions.bed")
    args = parser.parse_args()
    config = parseJSON(args.json)
    config.p = bool(args.p)
    config.c = bool(args.c)
    if config.p and config.c:
        abortWithMessage("Cannot declare -p and -c options together!")
    if config.p:
        config.outfile = 'Depth_of_Coverage-p.xlsx'
    elif config.c:
        config.outfile = 'Depth_of_Coverage-c.xlsx'
    else:
        config.outfile = 'Depth_of_Coverage.xlsx'
    if os.path.exists(config.outputDir) and config.outfile in [str(x) for x in os.listdir(config.outputDir)]:
        abortWithMessage("The directory {0} already exists and contains depth gauge output. Will not overwrite.".format(config.outputDir))
    elif not os.path.exists(config.outputDir):
        try:
            os.makedirs(config.outputDir)
        except OSError:
            abortWithMessage("Error when creating output directory {0}".format(config.outputDir))

    # check that all specified variant files exist
    for fn in config.inputFiles:
        if not os.path.exists(fn):
            abortWithMessage("Could not find variant file: {0}".format(fn))
    config.regions = parseAndValidateRegions(args.regions, config) 
    covDF = GaugeDepth(config)
    printOutput(config, covDF)
    # pretty print newline before exit
    print()


if __name__ == "__main__":
    if sys.hexversion < 0x02070000:
        raise RuntimeWarning("mucor should be run on python 2.7.0 or greater.")
    main()
