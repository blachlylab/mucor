#!/usr/bin/env python
#-*- coding: utf8
#
# mucor_config.py
#
# Karl Kroll
# James S Blachly

import os
import argparse
from collections import defaultdict
import csv
import json
import pdb
import codecs

# mucor modules
from config import Config

def abortWithMessage(message, help = False):
    print("*** FATAL ERROR: " + message + " ***")
    exit(2)

def throwWarning(message, help = False):
    print("*** WARNING: " + message + " ***")

def DetectDataType(fn):
    '''
    Pass a variant calls file, including path, and it will return the name of the tool that created the file. 
    Only works with the expected formats, outlined below.
    '''

    varFile = open(fn, 'r')
    varReader = csv.reader(varFile, delimiter='\t')
    row = varReader.next()
    while str(row).split("'")[1][0:2] == '##':
        if str('Torrent Unified Variant Caller') in str(row): 
            return "IonTorrent"
            break
        elif str('MiSeq') in str(row):
            return "MiSeq"
            break
        elif str('SomaticIndelDetector') in str(row):
            return "SomaticIndelDetector"
            break
        elif str('MuTect') in str(row):
            return "Mutect"
            break
        elif str('muTector') in str(row):
            return "muTector"
            break
        elif str('samtools') in str(row):
            return "Samtools"
            break
        elif str('source=VarScan') in str(row):
            return "VarScan"
            break
        elif str('ID=HaplotypeCaller') in str(row):
            return "HaplotypeCaller"
            break
        elif str('freeBayes') in str(row):
            return "FreeBayes"
            break
        elif str('source=GATK') in str(row):
            return "GenericGATK"
            break
        row = varReader.next()
    return "Unknown"

def DetectSnpEffStatus(fn):
    '''
    Pass a variant calls file, including path, and return True or False according to whether SnpEff has been run on this file.
    '''
    varFile = open(fn, 'r')
    varReader = csv.reader(varFile, delimiter='\t')
    try:
        row = varReader.next()
    except StopIteration:
        abortWithMessage("Empty File: " + str(fn))
    while str(row).split("'")[1][0:2] == '##':
        if str('SnpEff') in str(row): 
            return True
            break
        row = varReader.next()
    return False

def getJSONDict(args, proj_dir):
    '''
    Create the JSON dictionary. JSON dict keys are identical to the Config class in config.py. 
    Converts the user input args from argparse into a dictionary
    '''
    json_dict = defaultdict()
    json_dict['outputDir'] = str(args.output_directory).split('/')[-1]
    json_dict['gff'] = str(args.gff)
    json_dict['union'] = bool(args.union)
    json_dict['fast'] = args.archive_directory
    json_dict['feature'] = str(args.featuretype)
    json_dict['samples'] = list(dict())

    # VCF filters
    outFilters = set(["PASS"]) # By default, all mutations marked as PASS are permitted
    for i in str(args.vcf_filters).split(','):
        if i:
            outFilters.add(i) # Add all user-defined vcf filters. VCFs with empty 'FILTER' columns must have '.' supplied here, or all mutations will be filtered out
    json_dict['filters'] = [x for x in outFilters] # finally, convert set to list

    # Genomic regions
    json_dict['regions'] = []
    if args.regions: # user has defined regions to focus on
        for i in str(args.regions).split(','):
            if str(i.split('.')[-1]).lower() == "bed":
                # this is a bed file of regions
                # it will be parsed in the main python script
                if os.path.isfile(str(i)):
                    json_dict['regions'].append(os.path.expanduser(str(i)))
                else:
                    abortWithMessage("BED file {0} cannot be found".format(str(i)))
            elif str(str(i).split(':')[0]).startswith('chr') and str(str(i).split(':')[0]) != "chr":
                # this looks like a 'chromosome:start-stop' formatted region
                try:                                    # does the input region have valid start and ends?
                    int(str(str(i).split(':')[1])[0])
                    int(str(str(i).split(':')[1])[1])
                    json_dict['regions'].append(i)
                except ValueError:                      # start and/or end are invalid
                    abortWithMessage("Region {0} is not valid. Follow standard convention, Ex: chr1:100-300".format(str(i)))
                except IndexError:                      # only chromosome was defined. this is permitted (whole chromosome region)
                    json_dict['regions'].append(i)
            else:
                abortWithMessage("Region {0} is not a bed file or valid region.".format(i))

    # Variant databases 
    json_dict['databases'] = defaultdict(str)
    for i in args.databases:
        source_name = i.split(':')[0]
        source_path = i.split(':')[1]
        if bool(json_dict['databases'][source_name]):
            abortWithMessage("\ndatabases cannot share the same name.\nrename {0}:{1}\n".format(str(source_name), str(source_path)))
        json_dict['databases'][source_name] = os.path.expanduser(source_path)

    # Output formats
    json_dict['outputFormats'] = []
    for i in str(args.output_type).split(','):
        json_dict['outputFormats'].append(str(i).lower())

    # Samples and associated variant files
    for id in open(args.samples):
        sid = id.strip()
        if str(sid) == "":
            continue
        else:
            # sample name is not blank line
            # walk over the input project directory (or CWD if none defined) to find all files with sample ID in the name
            tmpSampleDict = defaultdict()
            tmpSampleDict['id'] = str(sid)
            tmpSampleDict['files'] = list()
            for root, dirs, files in os.walk(proj_dir):
                for i in files:
                    if str(sid) in str(i): # be careful with sample names here. "U-23" will catch "U-238" etc. Occasional cases can be resolved by manually editing the JSON config file
                        full_path = os.path.join(root, i)
                        source = DetectDataType(full_path)
                        if source == "Unknown":
                            print("Cannot parse file from an unsupported variant caller. \nPlease use supported variant software, or compose an input module for inputs.py")
                        if str(i).split('.')[-1] == str("vcf"):
                            tmpSampleDict['files'].append({'type':'vcf', 'path':str(full_path), 'snpeff':DetectSnpEffStatus(full_path), 'source':source} )
                        elif str(i).split('.')[-1] == str("out"):
                            tmpSampleDict['files'].append({'type':'mutect', 'path':str(full_path), 'source':source} )
                    # Not sure if these still work 
                        elif str(i).split('.')[-1].lower() == str("maf"):
                            tmpSampleDict['files'].append({'type':'maf', 'path':str(full_path), 'source':source} )
                        elif str(i).split('.')[-1].lower() == str("gvf"):
                            tmpSampleDict['files'].append({'type':'gvf', 'path':str(full_path), 'source':source} )
                        else:
                            # If not a VCF, MAF, GVF, or Mutect .out type, ignore it. Uncomment the following line to see the names of files that are being ignored
                            # print("Found an unsupported file type " + str(full_path) + " for sample " + str(sid))
                            pass
                            
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

def main():
    print
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gff", required=True, help="Annotation GFF/GTF for feature binning")
    parser.add_argument("-db", "--databases", default=[], action='append', help="Colon delimited name and path to variant database in bgzipped VCF format. Can be declared more than once. Ex: -db name1:/full/user/path/name1.vcf.gz ")
    parser.add_argument("-s", "--samples", required=True, help="Text file containing sample names")
    parser.add_argument("-d", "--project_directory", required=False, help="Working/project directory, in which to find output")
    parser.add_argument("-f", "--featuretype", required=True, help="Feature type into which to bin [gene]")
    parser.add_argument("-vcff", "--vcf_filters", default='', help="Comma separated list of VCF filters to allow")
    parser.add_argument("-a", "--archive_directory", default=False, help="Specify directory in which to read/write archived annotations. Undefined will circumvent annotation archive features.")
    parser.add_argument("-r", "--regions", default='', help="Comma separated list of bed regions OR bed files")
    # does union even work? should we delete the option entirely?
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
    parser.add_argument("-jco", "--json_config_output", required=True, help="Name of JSON configuration file")   
    parser.add_argument("-outd", "--output_directory", required=True, help="Name of Mucor output directory")
    parser.add_argument("-outt", "--output_type", default="default", help="Comma separated list of disired output types: xls, bed, long ")
    args = parser.parse_args()

    # Verify that required files and directoreis exist and/or are writable
    if not os.path.exists(args.gff):
        abortWithMessage("Could not find GFF file {0}".format(args.gff))
    if args.databases:
        for db in args.databases:
            try:
                if not os.path.exists(os.path.expanduser(db.split(':')[1])):
                    abortWithMessage("Could not find SNV DB file {0}".format(db.split(':')[1]))
            except IndexError: # user did not give a name and path, separated by colon
                abortWithMessage("Cannot process {0}\n\tDatabase input must be colon delimited as, 'database_name:database_path'".format(db))
    if not args.project_directory or not os.path.exists(args.project_directory):
        print("Project directory not found; using CWD")
        proj_dir = os.getcwd()
    else:
        proj_dir = args.project_directory
    if os.path.exists(args.json_config_output):
        abortWithMessage("JSON config file {0} already exists.".format(args.json_config_output))
    if os.path.exists(args.output_directory) and os.listdir(args.output_directory):
        abortWithMessage("The directory {0} already exists and contains output. Will not overwrite.".format(args.output_directory))
    elif not os.path.exists(args.output_directory):
        try:
            os.makedirs(args.output_directory)
        except:
            abortWithMessage("Error when creating output directory {0}".format(args.output_directory))

    # Construct JSON dictionary and dump it to the config file
    json_dict = getJSONDict(args, proj_dir)
    if inconsistentFilesPerSample(json_dict):
        throwWarning("Inconsistent number of files per sample")
        print("File #\tSample ID")
        for sample in json_dict['samples']:
            print( str(len(sample['files'])) + "\t" + str(sample['id']) )  
    output_file = codecs.open(args.json_config_output, "w", encoding="utf-8")
    json.dump(json_dict, output_file, sort_keys=True, indent=4, ensure_ascii=True)
    print

if __name__ == "__main__":
	main()
