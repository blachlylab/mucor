from __future__ import print_function
from pdb import set_trace as stop
import sys
import pandas as pd 
import os

from output import Writer
from config import Config

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

def main(fn):
    config = Config()
    config.outputDirName = './'
    print("parsing {0} ... ".format(fn))
    df = pd.io.excel.read_excel(fn)
    print("splitting ... ".format(fn))
    fixed = df.groupby(['chr','pos','ref','alt','feature','sample']).apply(unsplit)
    fixed.reset_index(drop=True,inplace=True)
    fixed.vf = fixed.vf.astype(float)
    ow = Writer()
    ow.write(fixed,'eucdist',config.outputDirName,config)

if __name__ == "__main__":
    if os.path.exists('vaf_euclidean_distance.xlsx'):
        print("will not overwrite ./feature_and_mutation_by_sample_vaf.xlsx\nmove it or rename it.")
        exit()

    main(sys.argv[1])