#!/usr/bin/env python
#-*- coding: utf8
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

#
# mucor_config.py
#
# Karl Kroll
# James S Blachly

import os
import sys
import re
import argparse
from collections import defaultdict
import csv
import json
import codecs
try:
    import HTSeq
except ImportError as err:
    print(str(err) + ". This module is required.")
    sys.exit()

# mucor modules
from inputs import Parser 
from config import Config
import output

def abortWithMessage(message, help = False):
    print("*** FATAL ERROR: " + message + " ***")
    exit(2)

def throwWarning(message, help = False):
    print("\n*** WARNING: " + message + " ***")

def DetectMalformedColumns(fn):
    '''
    Pass a variant calls file, including path, and it returns True if the column values are not consistent with VCF format
    Update: May not be required, now that HTSeq.VCF_Reader is used to parse
    '''
    varFile = open(fn, 'r')
    varReader = csv.reader(varFile, delimiter='\t')
    try:
        row = varReader.next()
    except StopIteration:
        throwWarning("Empty file: {0}".format(fn))
        return True
    while row[0].startswith('##'):
        try:
            row = varReader.next()
        except StopIteration:
            return True
    if row[:9] != ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
        return True
    return False


def DetectDataType(fn):
    '''
    Pass a variant calls file, including path, and it will return the name of the tool that created the file. 
    Only works with the expected formats, outlined below.
    Leverages the fact that most variant calling softwares output the name of the tool in every VCF header.
    '''
    unknown_data_type_name = "Unknown" # what to return when the variant caller cannot be inferred from the header
    try:
        varReader = HTSeq.VCF_Reader(str(fn))
        varReader.parse_meta() # extract list of sample names in VCF header
    except IOError:
        throwWarning("File cannot be read, or does not exist: {0}".format(fn))
        return unknown_data_type_name
    except KeyError:
        throwWarning("VCF file cannot be read. Is it malformed? {0}".format(fn))
        return unknown_data_type_name
    # HTSeq.VCF_Reader will accept and attempt to parse existing files that aren't VCF
    #   Parseing completes without warning, but the problem manifests as missing metadata (as if no VCF header)
    meta_info = str(varReader.meta_info()) 
    if not meta_info:
        if os.path.splitext(varReader.fos)[1] == ".vcf" or varReader.fos.endswith(".vcf.gz"):
            throwWarning("VCF file has no metadata. Is this a VCF with a header? {0}".format(fn))
            return unknown_data_type_name
        elif os.path.splitext(varReader.fos)[1] == ".out":
            return "muTector"
    if str('Torrent Unified Variant Caller') in meta_info: 
        return "IonTorrent"
    elif str('MiSeq') in meta_info:
        return "MiSeq"
    elif str('SomaticIndelDetector') in meta_info:
        return "SomaticIndelDetector"
    elif str('MuTect') in meta_info:
        return "Mutect"
    elif str('samtools') in meta_info:
        return "Samtools"
    elif str('VarScan') in meta_info:
        return "VarScan"
    elif str('HaplotypeCaller') in meta_info:
        return "HaplotypeCaller"
    elif str('freeBayes') in meta_info:
        return "FreeBayes"
    elif str('GATK') in meta_info: # this check is last, in case a more specific GATK module was detected above
        return "GenericGATK"
    else:
        return unknown_data_type_name

def DetectSnpEffStatus(fn):
    '''
    Pass a variant calls file, including path, and return True or False according to whether SnpEff has been run on this file.
    '''
    try:
        varReader = HTSeq.VCF_Reader(str(fn))
        varReader.parse_meta() # extract list of sample names in VCF header
    except IOError:
        throwWarning("File cannot be read, or does not exist: {0}".format(fn))
    except KeyError:
        throwWarning("VCF file cannot be read. Is it malformed? {0}".format(fn))
    metadata = varReader.metadata
    if str("SnpEff") in str("metadata"):
        return True
    else:
        return False

def blankJSONDict():
    '''
    Create a blank, valid JSON dictionary to show an example config file. JSON dict keys are identical to the Config class in config.py. 
    '''
    json_dict = defaultdict()
    json_dict['outputDir'] = str("./outputDir_path")
    json_dict['gff'] = str("~/ref/gff_path.gff")
    json_dict['union'] = bool(True)
    json_dict['fast'] = str("~/ref/fastDir_path/") # This will be boolean "False" by default, or a str() if declared
    json_dict['feature'] = str("gene_name")
    json_dict['samples'] = list(dict())

    # VCF filters
    json_dict['filters'] = ["PASS", "lowDP"] 

    # Genomic regions
    json_dict['regions'] = []

    # Variant databases 
    json_dict['databases'] = {"dbName":"~/ref/db_Location.vcf.gz"}

    # Output formats
    json_dict['outputFormats'] = []

    # Samples and associated variant files
    for sid in ['sample1','sample2']:
        # walk over the input project directory (or CWD if none defined) to find all files with sample ID in the name
        tmpSampleDict = defaultdict()
        tmpSampleDict['id'] = str(sid)
        tmpSampleDict['files'] = list()
        for root, dirs, files in ["~/project", "dir1", "sample1-a.vcf"], ["~/project", "dir1", "sample1.bam"], ["~/project", "dir2", "sample1-b.vcf"], ["~/project", "dir3", "sample2.vcf"], ["~/project", "dir3", "sample2.bam"]:
            if str(sid) in str(files): # be careful with sample names here. "U-23" will catch "U-238" etc. Occasional cases can be resolved by manually editing the JSON config file
                full_path = os.path.join(root, dirs, files)
                source = "VarScan"
                if os.path.splitext(files)[1] == str(".vcf") or files.endswith(".vcf.gz"):
                    tmpSampleDict['files'].append({'type':'vcf', 'path':str(full_path), 'snpeff':False, 'source':source} )
                elif os.path.splitext(files)[1] == str(".out"):
                    tmpSampleDict['files'].append({'type':'mutect', 'path':str(full_path), 'source':source} )
                elif os.path.splitext(files)[1] == str(".bam"):
                    try:
                        tmpSampleDict['files'].append({'type':'bam', 'path':str(full_path)} )
                    except KeyError:
                        tmpSampleDict['files'] = list() 
                        tmpSampleDict['files'].append({'type':'bam', 'path':str(full_path)})
            # Not sure if these still work 
                elif os.path.splitext(files)[1].lower() == str(".maf"):
                    tmpSampleDict['files'].append({'type':'maf', 'path':str(full_path), 'source':source} )
                elif os.path.splitext(files)[1].lower() == str(".gvf"):
                    tmpSampleDict['files'].append({'type':'gvf', 'path':str(full_path), 'source':source} )
                else:
                    # If not a VCF, MAF, GVF, or Mutect .out type, ignore it. Uncomment the following line to see the names of files that are being ignored
                    # print("Found an unsupported file type " + str(full_path) + " for sample " + str(sid))
                    pass
                            
        json_dict['samples'].append(tmpSampleDict)
    return json_dict

def processFile(full_path, tmpSampleDict):
    source = DetectDataType(full_path)
    if source == "Unknown" and (os.path.splitext(full_path)[1] == str(".vcf") or full_path.endswith(".vcf.gz") ):
        throwWarning(full_path)
        print("Cannot parse file from an unsupported or unknown variant caller. \nPlease use supported variant software, or compose an input module compatible with inputs.py")
        print("Options include: " + str(Parser().supported_formats.keys()).replace("'",""))
    elif os.path.splitext(full_path)[1] == str(".vcf") or full_path.endswith(".vcf.gz"):
        '''
        # may be handled by HTSeq.VCF_Reader 
        if DetectMalformedColumns(full_path):
            throwWarning(full_path)
            print("Malformed or empty VCF file: {0}".format(full_path))
        '''
        tmpSampleDict['files'].append({'type':'vcf', 'path':str(full_path), 'snpeff':DetectSnpEffStatus(full_path), 'source':source} )
    elif os.path.splitext(full_path)[1] == str(".out"):
        tmpSampleDict['files'].append({'type':'mutect', 'path':str(full_path), 'source':source} )
    elif os.path.splitext(full_path)[1] == str(".bam"):
        try:
            tmpSampleDict['files'].append({'type':'bam','path':str(full_path)} )
            # tmpSampleDict['files'].append({'type':'bam', 'path':str(full_path)} )
        except KeyError:
            tmpSampleDict['files'] = list()
            tmpSampleDict['files'].append({'type':'bam','path':str(full_path)} )
# Not sure if these still work 
    elif os.path.splitext(full_path)[1].lower() == str(".maf"):
        tmpSampleDict['files'].append({'type':'maf', 'path':str(full_path), 'source':source} )
    elif os.path.splitext(full_path)[1].lower() == str(".gvf"):
        tmpSampleDict['files'].append({'type':'gvf', 'path':str(full_path), 'source':source} )
    else:
        # If not a VCF, MAF, GVF, or Mutect .out type, ignore it. Uncomment the following line to see the names of files that are being ignored
        # print("Found an unsupported file type " + str(full_path) + " for sample " + str(sid))
        pass 
    return tmpSampleDict           

def parseAndValidateRegions(regions, json_config=False): 
    '''
    Read in a list of bed file locations and/or string-type chromosome ranges, uniquify it, and return a list of tuples in the form of (chrom, start, stop)
    Note: json_config is an optional argument. It is not used by mucor_config, but will be used by depth_gauge. This is because depth gauge can take regions via the passed json config and from the command line simultaneously. 
    '''
    if json_config:
        unique_region_list = [ x for x in set(regions.split(',')) | set(json_config.regions) if bool(x) ]
        verbose = True 
    else:
        unique_region_list = [ x for x in set(regions.split(',')) if bool(x) ]
        verbose = False # do not print running progress updates if executed from this mucor_config script
    regions_list = []
    if verbose:
        print("== Parsing and Validating Regions ==")
    for i in unique_region_list:
        # the parse
        if os.path.splitext(i)[1].lower() == ".bed":
            # this is a bed file of regions
            # it will be parsed in the main python script
            if os.path.isfile(str(i)):
                regions_list.append(os.path.expanduser(str(i)))
            else:
                throwWarning("BED file {0} cannot be found".format(str(i)))
                continue
        elif str(str(i).split(':')[0]).startswith('chr') and str(str(i).split(':')[0]) != "chr": # these two statements could be combined with a regular expression
            # this looks like a 'chromosome:start-stop' formatted region
            chrom = str(i).split(':')[0]
            start = 0
            end = 0
            try:                                    # does the input region have valid start and ends?
                start = int(re.split('-|\.\.', i.split(':')[1])[0]) # support chr:start-stop and chr:start..stop style location notation
                end   = int(re.split('-|\.\.', i.split(':')[1])[1])
                regions_list.append((str(chrom), start, end))
            except ValueError:                      # start and/or end are invalid, i.e. could not be converted to int
                throwWarning("Region {0} is not valid. Follow standard convention, Ex: chr1:100-300".format(str(i)))
                continue
            except IndexError:                      # only chromosome was defined, or only chromosome:start . this is permitted (whole chromosome region, or single base pair)
                if start and not end:
                    # single point
                    end = start
                    regions_list.append((str(chrom), start, end))
                elif not start and not end:
                    # whole contig
                    regions_list.append((str(chrom), None, None))
        else:
            throwWarning("Region {0} is not a bed file or valid region.".format(i))
        if verbose:
            print("\t{0}".format(i))

    if not regions_list:
        abortWithMessage("Regions not set")
    return regions_list

def getJSONDict(args):
    '''
    Create the JSON dictionary. JSON dict keys are identical to the Config class in config.py. 
    Converts the user input args from argparse into a dictionary
    '''
    json_dict = defaultdict()
    json_dict['outputDir'] = os.path.abspath(os.path.expanduser(str(args['output_directory'])))
    json_dict['gff'] = str(args['gff'])
    json_dict['union'] = bool(args['union'])
    json_dict['fast'] = args['archive_directory']
    json_dict['feature'] = str(args['featuretype'])
    json_dict['samples'] = list(dict())

    # VCF filters
    outFilters = set(["PASS", "."]) # By default, all mutations marked as PASS and '.' are permitted
    for i in str(args['vcf_filters']).split(','):
        if i:
            outFilters.add(i) # Add all user-defined vcf filters. VCFs with empty 'FILTER' columns must have '.' supplied here, or all mutations will be filtered out
    json_dict['filters'] = [x for x in outFilters] # finally, convert set to list

    # Genomic regions
    if args['regions']: # user has defined regions to focus on
        json_dict['regions'] = parseAndValidateRegions(args['regions'])
    else:
        json_dict['regions'] = args['regions']

    # Variant databases 
    json_dict['databases'] = defaultdict(str)
    for i in args['databases']:
        source_name = i.split(':')[0]
        source_path = i.split(':')[1]
        if bool(json_dict['databases'][source_name]):
            abortWithMessage("\ndatabases cannot share the same name.\nrename {0}:{1}\n".format(str(source_name), str(source_path)))
        json_dict['databases'][source_name] = os.path.expanduser(source_path)

    # Output formats
    json_dict['outputFormats'] = []
    for i in str(args['output_type']).split(','):
        json_dict['outputFormats'].append(str(i)) #.lower() was here for some reason

    # Samples and associated variant files
    samples = [] 
    for id in open(args['samples']):
        sid = id.strip()
        if str(sid) == "":
            continue
        else:
            # sample name is not blank line
            samples.append(sid)
    for sid in samples:
        # walk over the input project directory (or CWD if none defined) to find all files with sample ID in the name
        tmpSampleDict = defaultdict()
        tmpSampleDict['id'] = str(sid)
        tmpSampleDict['files'] = list()
        # directory crawl searching for sample name in the file name/path
        if args['project_directory']:
            for root, dirs, files in os.walk(args['project_directory']):
                for i in files:
                    if re.search(r"\b" + re.escape(sid) + r"\b", i): 
                    # be careful with sample names here. "U-23" will catch "U-238" etc. Occasional cases can be resolved by manually editing the JSON config file
                        full_path = os.path.abspath(os.path.join(root, i))
                        tmpSampleDict = processFile(full_path, tmpSampleDict) 
            #json_dict['samples'].append(tmpSampleDict)
        # if files were passed via --inputs
        if args['inputs']:
            for fn in args['inputs']:
                full_path = os.path.abspath(fn)
                if os.path.splitext(fn)[1] == str(".bam"):
                    if re.search(r"\b" + re.escape(sid) + r"\b", fn):
                        tmpSampleDict = processFile(full_path, tmpSampleDict)
                varReader = HTSeq.VCF_Reader(str(fn))
                varReader.parse_meta() # extract list of sample names in VCF header
                if sid in set(varReader.sampleids): 
                    tmpSampleDict = processFile(full_path, tmpSampleDict)
        # uniquify the list of files associated with the sample, in case they were supplied directly via (-i) and found while directory crawling (-d)
        tmpSampleDict['files'] =  { v['path']:v for v in tmpSampleDict['files']}.values()
        json_dict['samples'].append(tmpSampleDict)
    return json_dict

def inconsistentFilesPerSample(json_dict):
    '''
    For each sample, count the number of files found. If any sample has an inconsistent number compared to the others, returns True
    Useful for identifying sample IDs missing an expected file, or sample IDs that erroneously picked up extra files
    '''
    numset = set()
    for sample in json_dict['samples']:
        numset.add(len(sample['files']))
    if len(numset) > 1:
        return True

def filesWithoutSamples(json_dict, inputs):
    '''
    Survey the list of input files (explicitly passed with -i) and return those that are not part of the json dictionary.
    This indicates that the user wanted to analyze them with Mucor, but they were not appropriately identified/associated with a sample ID.
    '''
    # make a list of unique, full file paths associated with each sample
    #   i.e. the files that will be analyzed in the main script
    #   WARNING: complicated, chained python list comprehension ahead! 
    identified_file_paths = set([ os.path.abspath(x['path']) for b in json_dict['samples'] for x in b['files'] ])
    # make a list of unique, full file paths that were input to the config script via (-i)
    input_file_paths =  set([ os.path.abspath(x) for x in inputs ])
    # subtract the sets to identify the files passed via (-i) that were NOT found to be associated with any input sample ID
    missing_file_paths = input_file_paths - identified_file_paths
    return missing_file_paths

class exampleJSON(argparse.Action):
    '''
    This argparse action will run the example JSON config function and exit the program. It is executed when the user gives "-ex || --example" as an argument
    '''
    def __init__(self,option_strings, dest=None, nargs=0, default=None, required=False, type=None, metavar=None, help="Print a valid, example JSON config file and exit."):
        super(exampleJSON, self).__init__(option_strings=option_strings, dest=dest, nargs=nargs, default=default, required=required, metavar=metavar, type=type, help=help)
    def __call__(self, parser, namespace, values, option_string=None):
        print("Printing example_config.json ... ")
        output_file = codecs.open('example_config.json', "w", encoding="utf-8")
        json.dump(blankJSONDict(), output_file, sort_keys=True, indent=4, ensure_ascii=True)
        print("Finished. Exiting ... ")
        sys.exit(0)

def main():
    print
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-ex", "--example", action=exampleJSON)
    parser.add_argument("-g", "--gff", required=True, help="Annotation GFF/GTF for feature binning")
    parser.add_argument("-f", "--featuretype", required=True, help="Feature type into which to bin. Gencode GTF example: gene_name, gene_id, transcript_name, transcript_id, etc. ")
    parser.add_argument("-db", "--database", dest='databases', metavar='<dbName:/path/database.vcf.gz>', default=[], action='append', help="Colon delimited name and path to variant database in bgzipped VCF format. Can be declared >= 0 times.")
    parser.add_argument("-s", "--samples", metavar='<sample_list.txt>', required=True, help="Text file containing sample names. One sample per line.")
    parser.add_argument("-d", "--project_directory", metavar='<dirname>', required=False, help="Working/project directory, in which to find input variant call files.")
    parser.add_argument("-vcff", "--vcf_filters", default='', help="Comma separated list of VCF filters to allow. Default: PASS") # the defualt value is applied later on in the getJSONDict function, not here.
    parser.add_argument("-a", "--archive_directory", default="", help="Specify directory in which to read/write archived annotations. Undeclared will prevent using the annotation archive features.")
    parser.add_argument("-r", "--regions", default=[], help="Comma separated list of bed regions and/or bed files by which to limit output. Ex: chr1:10230-10240,chr2,my_regions.bed")
    parser.add_argument("-i", "--inputs", nargs="+", help="Input files")
    parser.add_argument("-u", "--union", action="store_true", help="""
        Join all items with same ID for feature_type (specified by -f)
        into a single, continuous bin. For example, if you want intronic
        variants counted with a gene, use this option. 
        WARNING, this will lead to spurious results for features that are duplicated on the same contig.
        When feature names are identical, the bin will range from the beginning of the first instance to the end of the last, even if they are several megabases apart.
        Refer to the documentation for a resolution using 'detect_union_bin_errors.py'
        """)
    parser.add_argument("-jco", "--json_config_output", required=True, help="Name of JSON configuration output file")   
    parser.add_argument("-outd", "--output_directory", required=True, help="Name of directory in which to write Mucor output")
    parser.add_argument("-outt", "--output_type", default="default", help=str("Comma separated list of desired output types. Options include: " + str(output.Writer().supported_formats.keys()).replace("'","") + ". Default: counts,txt"))
    args = vars(parser.parse_args()) # convert the argparser.Namespace to a dictionary, so values can be overwritten 
    
    # Verify user inputs
    # Is there a project directory defined? If so, does it exist?
    if args['project_directory'] and not os.path.exists(args['project_directory']):
        throwWarning("Supplied project directory does not exist: {0}".format(args['project_directory']))
        args['project_directory'] = None
    # Do the intput files exist?
    if args['inputs']:
        for fn in args['inputs']:
            if not os.path.exists(fn):
                throwWarning("Supplied file does not exist: {0}".format(fn))
                # remove this file from inputs
                args['inputs'] = [x for x in args['inputs'] if str(x) != str(fn)]
    # Did the user supply at least 1 valid input file or project directory?
    if not args['project_directory'] and not args['inputs']:
        abortWithMessage("Must supply a valid input file(s) (-i) and/or project directory (-d)")
    # Does the gtf exist?
    if not os.path.exists(args['gff']):
        abortWithMessage("Could not find GFF file {0}".format(args['gff']))

    # Do the database VCFs exist, and are they named properly?
    if args['databases']:
        for db in args['databases']:
            try:
                if not os.path.exists(os.path.expanduser(db.split(':')[1])):
                    abortWithMessage("Could not find SNV DB file {0}".format(db.split(':')[1]))
            except IndexError: # user did not give a name and path, separated by colon
                abortWithMessage("Cannot process {0}\n\tDatabase input must be colon delimited as, 'database_name:database_path'".format(db))

    # Will this output overwrite an existing JSON config file?
    if os.path.exists(args['json_config_output']):
        abortWithMessage("JSON config file {0} already exists.".format(args['json_config_output']))

    # Does the given output directory exist and contain output already?
    if os.path.exists(args['output_directory']) and [x for x in os.listdir(args['output_directory']) if x in output.Writer().file_names.values() ]:
        abortWithMessage("The directory {0} already exists and contains output. Will not overwrite.".format(args['output_directory']))
    elif not os.path.exists(args['output_directory']):
        # If it does not exist, is the parent directory writable?
        if not os.access( os.path.split( os.path.abspath(args['output_directory']) )[0], os.W_OK):
            abortWithMessage("Output directory {0} not writable".format(args['output_directory']))

    # Construct JSON dictionary and dump it to the config file
    json_dict = getJSONDict(args)
    if args['inputs']:
        missing_file_paths = filesWithoutSamples(json_dict, args['inputs']) 
        if missing_file_paths:
            throwWarning("Some input files were not associated with any sample name")
            for fn in missing_file_paths:
                print(fn)
    if inconsistentFilesPerSample(json_dict):
        throwWarning("Inconsistent number of files per sample")
        print("File #\tSample ID")
        for sample in json_dict['samples']:
            print( str(len(sample['files'])) + "\t" + str(sample['id']) )  
    output_file = codecs.open(args['json_config_output'], "w", encoding="utf-8")
    json.dump(json_dict, output_file, sort_keys=True, indent=4, ensure_ascii=True)
    print

if __name__ == "__main__":
	main()
