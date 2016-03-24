from __future__ import print_function
from pdb import set_trace as stop
import sys
import pandas as pd 
import os
import argparse 
import uuid

from output import Writer
from config import Config
import output

def unsplit(grp):
    # initialize empty variables for later
    outGrp = pd.DataFrame()

    # define useful variables
    # the columns we will ultimately want
    outCols = ['chr', 'pos', 'ref', 'alt', 'feature', 'sample'] # ['chr', 'pos', 'ref', 'alt', 'vf', 'dp', 'feature', 'effect', 'fc', '6500_Exomes', '1000_Genomes', '1000_Genomes_VAF', 'dbSNP-137-Common', 'Cosmic', 'dbSNP-142-All', 'dbSNP-142-Common', 'dbSNP-142-ClinVar', 'count', 'freq', 'sample', 'source']
    # problematic columns that need to be fixed
    badCols = [('VarScan_vf', 'VarScan_dp'), ('MuTect_vf',  'MuTect_dp')]
    # a "blank" output dataframe that has all info for the row, but 0 in vf
    blankGrp = grp[outCols]
    # the sources listed in the source column
    sources = grp.source.values[0].split(',')

    # loop over mutect- and varscan-specific columns, extracting respective vf when found
    for badCol in badCols: 
        if grp[badCol[0]].values != "?":
            _loopGrp = blankGrp
            # assign vf 
            _loopGrp['vf'] = grp[badCol[0]]
            #append this result to the output df
            outGrp = outGrp.append(_loopGrp)
    return outGrp.reset_index()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="the curated, sample_variant_details excel format file.")
    parser.add_argument("-outd", "--output_directory", required=True, help="Name of directory in which to write Mucor output")
    parser.add_argument("-outt", "--output_type", required=True, help=str("Comma separated list of desired output types. Options include: " + str(output.Writer().supported_formats.keys()).replace("'","") + ". Default: counts,txt"))
    args = parser.parse_args()

    fn = args.input
    config = Config()
    config.outputFormats.append(args.output_type)
    config.outputDirName = args.output_directory
    print("parsing {0} ... ".format(fn))
    sampleMask = {}
    for sheetname in ["calls with mutation type"]:
        df = pd.io.excel.read_excel(fn, sheetname=sheetname)
        df['sample'] = df['source']
        df.drop_duplicates(inplace=True)
        for samp in df['sample'].unique():
            if samp not in sampleMask.keys():
                sampleMask[str(samp)] = str(uuid.uuid4())
            df.replace(samp, sampleMask[samp], inplace=True)
        print("splitting ... ".format(fn))
        fixed = df.groupby(['chr','pos','ref','alt','feature','sample','VarScan_dp']).apply(unsplit)
        fixed.reset_index(drop=True,inplace=True)
        fixed.vf = fixed.vf.astype(float)
        ow = Writer()
        ow.write(fixed,args.output_type,config.outputDirName,config)
        thisName = os.path.join(config.outputDirName, 'feature_and_mutation_by_sample_vaf.xlsx')
        newName  = os.path.join(config.outputDirName, str(sheetname) + '.feature_and_mutation_by_sample_vaf.xlsx')
        os.rename(thisName, newName)
    pd.DataFrame(sampleMask, index=[0]).transpose().to_csv('sample_key.csv',sep=',', header=False, index=True)

if __name__ == "__main__":
    '''
    if os.path.exists('vaf_euclidean_distance.xlsx'):
        print("will not overwrite ./feature_and_mutation_by_sample_vaf.xlsx\nmove it or rename it.")
        exit()
    '''

    main()