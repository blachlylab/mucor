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

cwd = os.getcwd()

def abortWithMessage(message, help = False):
    print("*** FATAL ERROR: " + message + " ***")
    exit(2)

def throwWarning(message, help = False):
    print("*** WARNING: " + message + " ***")

def str_to_bool(s):
    if str(s) == 'False':
        return False
    else:
        return True

def DetectDataType(fn):
    #if str(fn).split('.')[-1].strip().lower() == str("maf"):
        #return "MAF"
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
        elif str('MuTect') in str(row) or str('muTector') in str(row):
            return "Mutect"
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
        row = varReader.next()
    return "Unknown"

def DetectSnpEffStatus(fn):
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

### TODO: snpEff trump detection? Ex: prefer a snp annotated version of a sample's output files

def filterFileList(json_dict_row):
    annotated_vcf = 3
    unannotated_vcf = 2
    mutect_output = 1
    status = defaultdict(dict)
    sid = json_dict_row['id']
    sources = []
    for i in json_dict_row['files']:
        if str(i['source']) not in sources:
            sources.append(str(i['source']))

    for i in json_dict_row['files']:
        if i['type'] == str('vcf') :
            if i['snpeff'] == bool(True) and annotated_vcf > status[i['source']][sid]:
                status[i['source']][sid] = annotated_vcf
            elif i['snpeff'] == bool(False) and unannotated_vcf > status[i['source']][sid]:
                status[i['source']][sid] = unannotated_vcf
        if i['type'] == str('mutect') and mutect_output > status[i['source']][sid]:
            mstatus[i['source']][sid] = mutect_output
    return status

def thing(args, proj_dir):
# TO DO : Pull these dictionary values from the 'Config' class variables for consistency
    json_dict = defaultdict()
    json_dict['outputDir'] = str(args.output_directory).split('/')[-1]
    json_dict['gff'] = str(args.gff)
    json_dict['union'] = bool(args.union)
    json_dict['fast'] = bool(args.no_archive)
    json_dict['feature'] = str(args.featuretype)
    json_dict['samples'] = list(dict())
    outFilters = set(["PASS"])
    for i in str(args.vcf_filters).split(','):
        if i:
            outFilters.add(i)
    json_dict['filters'] = [x for x in outFilters]
    json_dict['regions'] = []
    if args.regions:
        for i in str(args.regions).split(','):
            ## TO-DO: break these file/region checks out into separate, "isProperRegion" function(s)
            if str(i.split('.')[-1]).lower() == "bed":
                if os.path.isfile(str(i)):
                    json_dict['regions'].append(os.path.expanduser(str(i)))
                else:
                    abortWithMessage("BED file {0} cannot be found".format(str(i)))
            elif str(str(i).split(':')[0]).startswith('chr') and str(str(i).split(':')[0]) != "chr":
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
    json_dict['database'] = []
    for i in str(args.database).split(','):
        json_dict['database'].append(os.path.expanduser(i))
    json_dict['outputFormats'] = []
    for i in str(args.output_type).split(','):
        json_dict['outputFormats'].append(str(i).lower())
    for id in open(args.samples):
        sid = id.strip()
        if str(sid) == "":
            continue
        else:
            something = defaultdict()
            something['id'] = str(sid)
            something['files'] = list()
            for root, dirs, files in os.walk(proj_dir):
                for i in files:
                    if str(sid) in str(i): # be careful with sample names here. "U-24" will catch "U-240" etc. 
                        full_path = os.path.join(root, i)
                        if str(i).split('.')[-1] == str("vcf"):
                            something['files'].append({'type':'vcf', 'path':str(full_path), 'snpeff':DetectSnpEffStatus(full_path), 'source':DetectDataType(full_path)} )
                        elif str(i).split('.')[-1] == str("out"):
                            something['files'].append({'type':'mutect', 'path':str(full_path), 'source':DetectDataType(full_path)} )
                        elif str(i).split('.')[-1].lower() == str("maf"):
                            something['files'].append({'type':'maf', 'path':str(full_path), 'source':DetectDataType(full_path)} )
                        elif str(i).split('.')[-1].lower() == str("gvf"):
                            something['files'].append({'type':'gvf', 'path':str(full_path), 'source':DetectDataType(full_path)} )
                        else:
                            print("Found an unsupported file type " + str(full_path) + " for sample " + str(sid))
        json_dict['samples'].append(something)
    return json_dict

def inconsistentFilesPerSample(json_dict):
    numset = set()
    for sample in json_dict['samples']:
        numset.add(len(sample['files']))
    if len(numset) > 1:
        return True

def main():
    print
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gff", required=True, help="Annotation GFF/GTF for feature binning")
    parser.add_argument("-db", "--database", default=[], help="Comma separated list of known SNV databases in VCF format")
    parser.add_argument("-s", "--samples", required=True, help="Text file containing sample names")
    parser.add_argument("-d", "--project_directory", required=False, help="Working/project directory, in which to find output")
    parser.add_argument("-f", "--featuretype", required=True, help="Feature type into which to bin [gene]")
    parser.add_argument("-vcff", "--vcf_filters", default='', help="Comma separated list of VCF filters to allow")
    parser.add_argument("-n", "--no_archive", action="store_false", default=True, help="prevent quick load of annotation files")
    parser.add_argument("-r", "--regions", default='', help="Comma separated list of bed regions OR bed files")
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
    if not os.path.exists(args.gff):
        abortWithMessage("Could not find GFF file {0}".format(args.gff))
    if args.database:
        for db in str(args.database).split(','):
            if not os.path.exists(os.path.expanduser(db)):
                abortWithMessage("Could not find SNV DB file {0}".format(db))
    if not args.project_directory or not os.path.exists(args.project_directory):
        print("Project directory not found; using CWD")
        proj_dir = cwd
    else:
        proj_dir = args.project_directory
    json_dict = thing(args, proj_dir)
    if inconsistentFilesPerSample(json_dict):
        throwWarning("Inconsistent number of files per sample")
        print("File #\tSample ID")
        for sample in json_dict['samples']:
            print( str(len(sample['files'])) + "\t" + str(sample['id']) )
    if os.path.exists(args.json_config_output):
        abortWithMessage("JSON config file {0} already exists.".format(args.json_config_output))
    if os.path.exists(args.output_directory) and os.listdir(args.output_directory):
        abortWithMessage("The directory {0} already exists and contains output. Will not overwrite.".format(args.output_directory))
    elif not os.path.exists(args.output_directory):
        try:
            os.makedirs(args.output_directory)
        except:
            abortWithMessage("Error when creating output directory {0}".format(args.output_directory))
    output_file = codecs.open(args.json_config_output, "w", encoding="utf-8")
    json.dump(json_dict, output_file, sort_keys=True, indent=4, ensure_ascii=True)
    print

if __name__ == "__main__":
	main()
