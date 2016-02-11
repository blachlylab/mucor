import os
import pandas as pd
from config import Config
import numpy as np
import pdb
import time
from info import Info
import pdb


def VAFEuclideanDistance(self):
    '''
    Print n by n matrix, where n is the number of samples.
    Table values are the hyper-dimensional euclidean distance between samples,  
        where axes are "mutation" and values are "VAF".

    Output: vaf_euclidean_distance.xlsx
    '''
    self.attempted_formats.append('eucdist')
    outputDirName = self.outputDirName
    outputFileName = self.file_names['eucdist']
    varDF = self.data
    ofVAFEuclideanDistance = pd.ExcelWriter(str(outputDirName) + "/" + outputFileName)
    
    outDF = pd.DataFrame(varDF, columns=['feature','chr','pos','ref','alt','sample', 'vf'])
    outDF = pd.pivot_table(outDF, values='vf', index=['feature','chr','pos','ref','alt'], columns='sample')
    outDF.fillna(0,inplace=True)
    samps = set(outDF.columns)
    eucDF = pd.DataFrame(columns=samps, index=samps)
    while samps:
        sample = samps.pop()
        for other in samps:
            eucDist = np.sqrt((outDF[sample].subtract(outDF[other])**2).sum())
            eucDF.loc[sample, other] = eucDist
            eucDF.loc[other, sample] = eucDist
            # this copied, inverted-indeces line makes the output a symetric square
            # remove one of these to make output in the form of a non-redundant triangle
    eucDF.to_excel(ofVAFEuclideanDistance, 'Mutation VAF Distance Matrix', na_rep="", index=True)
    try:
        ofVAFEuclideanDistance.save()
    except IOError:
        abortWithMessage("Error opening output file {0}/{1}".format(outputDirName, outputFileName))
    print("\t{0}: {1} rows".format(str(outputDirName) + "/" + outputFileName, len(outDF)))        
    return True 