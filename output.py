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

# output.py
#
# James S Blachly, MD
# 2014-10-04

import os
import pandas as pd
from config import Config
import numpy as np
import time
from info import Info
from pdb import set_trace as stop

from mucor import abortWithMessage
from mucor import throwWarning

class Writer(object):
    """Object that parses the mucor dataframe and can write output in several different formats"""
    
    def __init__(self):
        self.data = pd.DataFrame()
        self.config = Config()
        self.outputDirName = ''
        self.attempted_formats = []     # used to prevent output modules from being executed multiple times
        self.supported_formats = {  "default": self.Default,
                                    "counts": self.Counts,
                                    "txt": self.VariantDetails,
                                    "longtxt": self.LongVariantDetails,
                                    "xls": self.VariantDetails,
                                    "longxls": self.LongVariantDetails,
                                    "bed":self.VariantBed,
                                    "featXsamp": self.FeatureXSample,
                                    "mutXsamp": self.Feature_and_MutationXSample,
                                    "mutXsampVAF": self.Feature_and_MutationXSampleVAF,
                                    "vcf": self.VCF,
                                    "all": self.All,
                                    "runinfo": self.RunInfo }

        self.file_names        = {  "counts": "counts.txt",
                                    "txt": "variant_details.txt",
                                    "longtxt": "long_variant_details.txt",
                                    "xls": "variant_details.xlsx",
                                    "longxls": "long_variant_details.xlsx",
                                    "bed": "variant_locations.bed",
                                    "featXsamp": "feature_by_sample.xlsx",
                                    "mutXsamp": "feature_and_mutation_by_sample.xlsx",
                                    "mutXsampVAF": "feature_and_mutation_by_sample_vaf.xlsx",
                                    "vcf": "variant_locations.vcf",
                                    "runinfo": "run_info.txt" }
        
    def write(self, data, format, outputDirName, config):
        '''
        Write data in format to outputDirName
        data: a pandas dataframe
        format: a string: one of [default, vcf, gvf, xls]   # consider merging in report functionality
        outputDirName: a string
        '''
        
        # Parameter integrity/validity checks
        if (type(data) == pd.DataFrame):
            self.data = data   
        else:
            raise TypeError("data was not pandas.core.frame.DataFrame")
        
        if (format not in self.supported_formats.keys()):
            raise ValueError("Unsupported output format {0}\n\tSupported Formats: {1}".format(format, self.supported_formats.keys()))
        
        if not os.path.exists(outputDirName):
            raise ValueError("The specified output directory ({0}) does not exist.".format(outputDirName))
        else:
            self.outputDirName = outputDirName
        self.format = format
        self.config = config
        self.supported_formats[format]()
        
    def RunInfo(self):
        '''
        Print useful information about the run, including the version, time, and configuration used 

        '''
        self.attempted_formats.append('runinfo')
        outputDirName = self.outputDirName
        outputFileName = self.file_names['runinfo']
        config = self.config
        try: 
            ofRunInfo = open(outputDirName + "/" + outputFileName, 'w+')
        except IOError:
            abortWithMessage("Error opening output file {0}/{1}".format(outputDirName, outputFileName))
        # =========================
        # run_info.txt
        #
        ofRunInfo.write(Info.versionInfo + "\n")
        ofRunInfo.write("{0}\n\n".format(time.ctime() ) )
        ofRunInfo.write(str(config))
        
        '''
        TO DO: total runtime? or runtime breakdown per segment of the program?
        ofRunInfo.write("Variants Pre-filter: \n")
        ofRunInfo.write("        Post-filter: \n")
        '''
        ofRunInfo.close()
        return True

    def Counts(self):
        '''
        Print counts per feature

        Output: counts.txt
        '''
        self.attempted_formats.append('counts')
        outputDirName = self.outputDirName
        outputFileName = self.file_names['counts']
        varDF = self.data
        try:
            ofCounts = open(outputDirName + "/" + outputFileName, 'w+')
        except IOError:
            abortWithMessage("Error opening output file {0}/{1}".format(outputDirName, outputFileName))
        
        grouped = varDF.groupby('feature')

        numHits = grouped['sample'].count()
        numHits.name = 'Hits'
        numHits = pd.DataFrame(numHits)
        
        weightedHits = grouped['vf'].apply(np.sum)
        weightedHits.name = 'WeightedHits'
        weightedHits = pd.DataFrame(weightedHits)

        avgWeight = weightedHits.div(numHits.Hits, axis='index')
        avgWeight = avgWeight.rename(columns={'WeightedHits':'AverageWeight'}) # different kind of column re-name required for division result
        avgWeight = pd.DataFrame(avgWeight)

        uniqueHits = grouped['pos'].nunique()
        uniqueHits.name = 'UniqueHits'
        uniqueHits = pd.DataFrame(uniqueHits)

        numSamples = grouped['sample'].nunique()
        numSamples.name = 'NumSamples'
        numSamples = pd.DataFrame(numSamples)

        # merge the 5 1-column dataframes into a single, 5-column dataframe. (6 columns including the feature name column)
        out = numHits.join(weightedHits).join(avgWeight).join(uniqueHits).join(numSamples)
        # change capitalization for consistency 
        out.index.names = ['FeatureName']

        mySort(out).to_csv(ofCounts, sep='\t', na_rep='?', index=True)
        print("\t{0}: {1} rows".format(ofCounts.name, len(out)))        
        return True 

    def FeatureXSample(self):
        '''
        Print counts per feature per sample.
        Sample names populate the table header and feature names populate the first column. The count per sample per feature are the table values.

        Output: feature_by_sample.xlsx
        '''
        self.attempted_formats.append('featXsamp')
        outputDirName = self.outputDirName
        outputFileName = self.file_names['featXsamp']
        varDF = self.data
        ofFeatureXSample = pd.ExcelWriter(str(outputDirName) + "/" + outputFileName)
        
        groupedDF = pd.DataFrame(varDF.groupby(['feature','sample']).apply(len))
        outDF = groupedDF.stack().unstack(1)
        outDF.index = outDF.index.droplevel(1)
        outDF.to_excel(ofFeatureXSample, 'Feature by Sample', na_rep=0, index=True)
        try:
            ofFeatureXSample.save()
        except IOError:
            abortWithMessage("Error opening output file {0}/{1}".format(outputDirName, outputFileName))
        print("\t{0}: {1} rows".format(str(outputDirName) + "/" + outputFileName, len(outDF)))        
        return True 

    def Feature_and_MutationXSample(self):
        '''
        Print counts per mutation per sample.
        Sample names populate the table header and feature names populate the first column. Chromosome, position, ref, and alt populate the next columns.
        The table values are boolean: 1 for present mutation, 0 for missing mutation. 

        Output: feature_and_mutation_by_sample.xlsx
        '''
        self.attempted_formats.append('mutXsamp')
        outputDirName = self.outputDirName
        outputFileName = self.file_names['mutXsamp']
        varDF = self.data
        ofFeature_and_MutationXSample = pd.ExcelWriter(str(outputDirName) + "/" + outputFileName)
        
        groupedDF = pd.DataFrame(varDF.groupby(['feature','chr','pos','ref','alt','sample']).apply(len))
        outDF = groupedDF.stack().unstack(5)
        outDF.index = outDF.index.droplevel(5)
        outDF.to_excel(ofFeature_and_MutationXSample, 'Feature and Mutation by Sample', na_rep=0, index=True)
        try:
            ofFeature_and_MutationXSample.save()
        except IOError:
            abortWithMessage("Error opening output file {0}/{1}".format(outputDirName, outputFileName))
        print("\t{0}: {1} rows".format(str(outputDirName) + "/" + outputFileName, len(outDF)))        
        return True 

    def Feature_and_MutationXSampleVAF(self):
        '''
        Print counts per mutation per sample.
        Sample names populate the table header and feature names populate the first column. Chromosome, position, ref, and alt populate the next columns.
        The table values are mutation VAF. 

        Output: feature_and_mutation_by_sample_vaf.xlsx
        '''
        self.attempted_formats.append('mutXsampVAF')
        outputDirName = self.outputDirName
        outputFileName = self.file_names['mutXsampVAF']
        varDF = self.data
        ofFeature_and_MutationXSample = pd.ExcelWriter(str(outputDirName) + "/" + outputFileName)
        
        outDF = pd.DataFrame(varDF, columns=['feature','chr','pos','ref','alt','sample', 'vf'])
        outDF = pd.pivot_table(outDF, values='vf', index=['feature','chr','pos','ref','alt'], columns='sample')
        outDF.to_excel(ofFeature_and_MutationXSample, 'Feat. + Mutation by Sample VAF', na_rep=0, index=True)
        try:
            ofFeature_and_MutationXSample.save()
        except IOError:
            abortWithMessage("Error opening output file {0}/{1}".format(outputDirName, outputFileName))
        print("\t{0}: {1} rows".format(str(outputDirName) + "/" + outputFileName, len(outDF)))        
        return True 

    def VariantDetails(self):
        '''
        Print all information about each mutation,
        combining all mutations (irrespective of in how many samples they appear)
        into a single row
        Note: chrom, position, ref, alt, and feature are all required to uniquely identify a mutation 
              indels may have the same chr, pos, but different ref/alt

        Output: variant_details.txt, variant_details.xls
        Note: switching the pandas ExcelWriter file extension to xlsx instead of xls requires openpyxl
              If openpyxl is unavailable, the xlwt library can write xls 
        '''

        outputDirName = self.outputDirName
        
        varDF = self.data
        if 'txt' in self.config.outputFormats and 'txt' not in self.attempted_formats: # and not os.path.exists(outputDirName + "/variant_details.txt"):
            txt = bool(True)
            self.attempted_formats.append('txt')
        else:
            txt = bool(False) # User has not opted for this output [parenthetical 'if's] or this output has already been run [last 'if']
        if 'xls' in self.config.outputFormats and 'xls' not in self.attempted_formats: # and not os.path.exists(outputDirName + "/variant_details.xls"):
            xls = bool(True)
            self.attempted_formats.append('xls')
        else:
            xls = bool(False) # User has not opted for this output [parenthetical 'if's] or this output has already been run [last 'if']

        try:
            if txt:
                outputFileName = self.file_names['txt']
                ofVariantDetailsTXT = open(outputDirName + "/" + outputFileName, 'w+')
            if xls:
                outputFileName = self.file_names['xls']
                ofVariantDetailsXLS = pd.ExcelWriter(str(outputDirName) + '/' + outputFileName)
        except IOError:
            abortWithMessage("Error opening output file {0}/{1}".format(outputDirName, outputFileName))
        if txt or xls:
            # Group by (chr, pos, ref, alt, feature)
            grouped = varDF.groupby(['chr', 'pos', 'ref', 'alt', 'feature'])
            # apply collapsing function to each pandas group
            outSeries = grouped.apply(collapseVariantDetails)
            out = pd.DataFrame(outSeries.reset_index(drop=True), columns=varDF.columns)

            #outSeries = grouped.apply(collapseVariantDetails)
            #out = pd.DataFrame(outSeries.reset_index(drop=True), columns=varDF.columns)

        if txt:
            # print the new, collapsed dataframe to a file
            mySort(out, ['feature','pos']).to_csv(ofVariantDetailsTXT, sep='\t', na_rep='?', index=False)
            print("\t{0}: {1} rows".format(ofVariantDetailsTXT.name, len(out)))
        if xls:
            # print the new, collapsed dataframe to file a
            mySort(out, ['feature','pos']).to_excel(ofVariantDetailsXLS, 'Variant Details', na_rep='?', index=False)
            ofVariantDetailsXLS.save()
            print("\t{0}: {1} rows".format(str(outputDirName + '/' + outputFileName), len(out)))

        return True

    def LongVariantDetails(self):
        '''
        Similar to printVariantDetails above, but writes each instance
        of a mutation to a new row. 
        Each mutation is written once per source instead of combining
        reoccurring mutations in to 1 unique row.
        
        Output: long_variant_details.txt, long_variant_details.xls
        Note: switching the pandas ExcelWriter file extension to xlsx instead
        of xls requires openpyxl
        '''
        outputDirName = self.outputDirName
        
        varDF = self.data

        if 'longtxt' in self.config.outputFormats and 'longtxt' not in self.attempted_formats:
            # and not os.path.exists(outputDirName + "/long_variant_details.txt"):
            longtxt = bool(True)
            self.attempted_formats.append('longtxt')
        else:
            longtxt = bool(False) # User has not opted for this output [parenthetical 'if's] or this output has already been run [last 'if']
        if 'longxls' in self.config.outputFormats and 'longxls' not in self.attempted_formats: # and not os.path.exists(outputDirName + "/long_variant_details.xls"):
            longxls = bool(True)
            self.attempted_formats.append('longxls')
        else:
            longxls = bool(False) # User has not opted for this output [parenthetical 'if's] or this output has already been run [last 'if']

        try:
            if longtxt:
                outputFileName = self.file_names['longtxt']
                ofLongVariantDetailsTXT = open(outputDirName + "/" + outputFileName, 'w+')
                mySort(varDF, ['feature','pos']).to_csv(ofLongVariantDetailsTXT, sep='\t', na_rep='?', index=False)
                ofLongVariantDetailsTXT.close()
                print("\t{0}: {1} rows".format(str(outputDirName + '/' + outputFileName), len(varDF)))
                    
            if longxls:
                outputFileName = self.file_names['longxls']
                ofLongVariantDetailsXLS = pd.ExcelWriter(str(outputDirName) + '/' + outputFileName)
                mySort(varDF, ['feature','pos']).to_excel(ofLongVariantDetailsXLS, 'Long Variant Details', na_rep='?', index=False)
                ofLongVariantDetailsXLS.save()
                print("\t{0}: {1} rows".format(str(outputDirName + '/' + outputFileName), len(varDF)))
        except IOError:
            abortWithMessage("Error opening output file {0}/{1}".format(outputDirName, outputFileName))
        
        return True

    def VariantBed(self):
        '''
        Print bed file of the variant locations

        Output: variant_locations.bed
        '''
        self.attempted_formats.append('bed')
        outputDirName = self.outputDirName
        outputFileName = self.file_names['bed']
        varDF = self.data
        try:
            ofVariantBeds = open(outputDirName + "/" + outputFileName, 'w+')
        except IOError:
            abortWithMessage("Error opening output file {0}/{1}".format(outputDirName, outputFileName))
        
        grouped = varDF.groupby(['chr', 'pos', 'ref', 'alt'])
        outSeries = grouped.apply(collapseVariantBed)
        out = pd.DataFrame(outSeries.reset_index(drop=True))
        mySort(out, ['chr','start']).to_csv(ofVariantBeds, sep='\t', na_rep='?', index=False, header=False)
        ofVariantBeds.close()
        print("\t{0}: {1} rows".format(ofVariantBeds.name, len(out)))
        return True

    def VCF(self):
        '''
        Print vcf file of the variant locations, feature, depth, and variant frequency.

        Output: variant_locations.vcf
        '''
        self.attempted_formats.append('vcf')
        outputDirName = self.outputDirName
        outputFileName = self.file_names['vcf']
        varDF = self.data
        try:
            ofVariantVCF = open(outputDirName + "/" + outputFileName, 'w+')
        except IOError:
            abortWithMessage("Error opening output file {0}/{1}".format(outputDirName, outputFileName))
         
        #VCF header to stream
        ofVariantVCF.write("##fileformat=VCFv4.1\n")
        
        grouped = varDF.groupby(['chr', 'pos', 'ref', 'alt'])
        vcf_fields = ['chr','pos','id','ref','alt','qual','filter','INFO']  # the relevant fields for a vcf file
        outSeries = grouped.apply(collapseVCF)                           # make the INFO column for this set of samples
        outDF = pd.DataFrame(outSeries, columns=['INFO'])
        varDF = varDF.merge(outDF, left_on=['chr','pos','ref','alt'], right_index=True) # add the info column to the original dataframe
        out = varDF.reindex(columns=vcf_fields).fillna('.')                 # drop the columns that are unrelated to vcf format, while reordering columns into proper vcf order
        official_vcf_fields = ['#CHROM','POS','ID', 'REF','ALT', 'QUAL', 'FILTER','INFO']   
        out.columns = pd.Index(official_vcf_fields)                         # rename columns to comply with official vcf standards
        out.to_csv(ofVariantVCF, sep='\t', na_rep='?', index=False, header=True, sparsify=False)
        ofVariantVCF.close()
        print("\t{0}: {1} rows".format(ofVariantVCF.name, len(out)))
        return True


    def Default(self):
        '''
        Runs several, common output functions, including VariantDetails and Counts. 
        RunInfo is executed from the main printOutput function, regardless of the output types the user has selected
        
        Output: variant_details.txt
                counts.txt
        '''
        for format in ['runinfo', 'txt', 'counts']:
            self.supported_formats[format]()
        return True        

    def All(self):
        '''
        Runs all available output functions, based on those in the self.supported_formats dictionary 
        '''
        formats = self.supported_formats.keys()
        formats.remove('default') # the default output is redundant in this situation, since we are going to run all the functions anwyay.
        formats.remove('all')     # remove 'all' or initiate an infinite loop
        self.config.outputFormats = formats
        for format in formats:
            self.supported_formats[format]()
        return True


############################################
#     |     support functions      |       #
#     V                            V       #
# do not write output, but used by writers #
############################################

def collapseVariantDetails(group):
    '''
    Pandas operations to support the VariantDetails family of functions. 
    Collapses variant rows that share the same contig, position, ref allele, alt allele, and feature.
    Input: a pandas groupby object
    Output: a pandas dataframe object
    '''
    outvals = []
    columns = list(group.keys().values)
    for column in columns:
        # outstring = '' # string based
        outlist = []     # list based
        if column in ['vf', 'dp', 'sample', 'source']:  # the only columns that need to be concatenated 
                                                        # the others are uniquified and "always" yield 1 value
                                                                                    
            for i in group[column].values:
                #outstring += str(i) + ", "   # string based
                outlist.append(str(i))
            #outvals.append( outstring[:-2] )  # string based          # trim the extra ', ' off the end of outstring 
            outvals.append(", ".join(outlist))
        else:
            if len(group[column].unique()) == 1:
                outvals.append( group[column].unique()[0] ) # only 1 value; extract it since there's no need to preserve the array type
            else:
                #this should be pretty rare, such as inputting two files for the same sample that were annotated for functional consequence by two different tools
                outvals.append( "/".join(filter(bool, [x for x in group[column].unique() if x != '?' ] )) )  # functional consequence and effect should only be 1 value as well, but may not be if there is a mixture of annotated and unannotated vcf files. 
                                                                    # This sorts the effects and consequences, and concatenates the values that are not '?'s.
                                                                    # each of the appended values must be an array in order to transform dictionary to pandas data frame in the next line
    try:                                                           
        #outDF = pd.DataFrame( dict(zip(columns,outvals)), columns=columns )
        outD = pd.Series( dict(zip(columns,outvals)) )
    except ValueError:
        abortWithMessage("Could not collapse variant {0} {1} {2}/{3} into a single row. \n\t\tThis mutation may have inconsistent effect and/or functional consequence values in different VCF files.".format( str(group['chr'].unique()[0]), str(group['pos'].unique()[0]), str(group['ref'].unique()[0]), str(group['alt'].unique()[0]) ))
    return outD

def collapseVariantBed(group):
    '''
    Pandas operation to support the VariantBed function.
    Collapses variant rows that share the same contig, position, ref allele, and alt allele
    Input: a pandas groupby object
    Output: a pandas dataframe object, ready to be written to bed format
    '''
    bed_fields = ['chr','start','end','name']
    chrom = group['chr'].unique()[0]
    pos = group['pos'].unique()[0]
    start = pos - 1     # bed files are 0-based positions
    end = pos
    name = ";".join(set(group['feature'].values))   # mutations that fall in overlapping features will have a name column with multiple values delimited by semicolons
    return pd.Series( dict( zip(     \
        ['chr','start','end','name'],   \
        [chrom, start , end , name]     \
        ))) # , columns=bed_fields

def collapseVCF(group):
    '''
    Pandas operation to support the VCF function.
    Collapses variant rows that share the same contig, position, ref allele, and alt allele
    Input: a pandas groupby object
    Output: a pandas dataframe object, ready to be written to vcf format 

    TODO: make this output a (multi-)sample vcf, i.e. each sample in its own column
    '''
    info = '' # output will have sample-specific information in the 'info' column
    samples = group['sample'].unique()
    for sample in samples:
        info += str(sample) + ':'
        feat = "/".join(group[group['sample'] == sample]['feature'].unique())
        lines = group[ (group['sample'] == sample) ]
        # pandas 0.16.1 added pandas.DataFrame.sample() which breaks
        # a query of form df.sample ==
        #line = group.query('(sample == @sample) and (feature == @feat)', engine='python') 
        vf = str(lines['vf'].unique()[0])
        dp = str(lines['dp'].unique()[0])
        info += str(feat) + ",vf=" + vf + ",dp=" + dp + ";"
    return info 

def mySort(df, columns=None):
    '''
    A wrapper for pandas.DataFrame sort that is agnostic about which version of pandas is used 
    '''
    version = float(".".join(str(pd.__version__).split('.')[0:2]))
    if version >= float(0.17):
        if columns:
            #sort by columns
            return df.sort_values(by=columns)
        else:
            #sort by index
            return df.sort_index()
    else:
        return df.sort(columns=columns)
