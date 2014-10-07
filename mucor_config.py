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

def str_to_bool(s):
    if str(s) == 'False':
        return False
    else:
        return True

def DetectDataType(fn):
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
        row = varReader.next()
    return "Unknown"

def DetectSnpEffStatus(fn):
    varFile = open(fn, 'r')
    varReader = csv.reader(varFile, delimiter='\t')
    row = varReader.next()
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
    json_dict['output_dir'] = str(args.output_directory).split('/')[-1]
    json_dict['gtf'] = str(args.gff)
    json_dict['union'] = bool(args.union)
    json_dict['fast'] = bool(args.no_archive)
    json_dict['feature'] = str(args.featuretype)
    json_dict['samples'] = list(dict())
    json_dict['output_formats'] = ['default', 'long']   # future: xls, vcf, gvf
    outFilters = set(["PASS"])
    for i in str(args.vcf_filters).split(','):
        if i:
            outFilters.add(i)
    json_dict['filters'] = [x for x in outFilters]
    json_dict['database'] = []
    for i in str(args.database).split(','):
        json_dict['database'].append(os.path.expanduser(i))
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
                        else:
                            print("unsure of what to do with " + str(full_path))
        json_dict['samples'].append(something)
    return json_dict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gff", required=True, help="Annotation GFF/GTF for feature binning")
    parser.add_argument("-db", "--database", default=[], help="Comma separated list of known SNV databases in VCF format")
    parser.add_argument("-s", "--samples", required=True, help="Text file containing sample names")
    parser.add_argument("-d", "--project_directory", required=False, help="Project root directory, in which to find output")
    parser.add_argument("-f", "--featuretype", required=True, help="Feature type into which to bin [gene]")
    parser.add_argument("-vcff", "--vcf_filters", default='', help="Comma separated list of VCF filters to allow")
    parser.add_argument("-n", "--no_archive", action="store_false", default=True, help="prevent quick load of annotation files")
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
    parser.add_argument("-od", "--output_directory", required=True, help="Name of Mucor output directory")
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
    if os.path.exists(args.json_config_output):
        abortWithMessage("JSON config file {0} already exists.".format(args.json_config_output))
    output_file = codecs.open(args.json_config_output, "w", encoding="utf-8")
    json.dump(json_dict, output_file, sort_keys=True, indent=4, ensure_ascii=True)

if __name__ == "__main__":
	main()
