#
# output.py
#
# James S Blachly, MD
# 2014-10-04

import os
import pandas as pd
from config import Config
import numpy as np
import pdb
import time
from info import Info

class Writer(object):
    """Object that parses the mucor dataframe and can write output in several different formats"""
    
    def __init__(self):
        self.data = pd.DataFrame()
        self.config = Config()
        self.outputDirName = ''
        self.supported_formats = {  "default": self.Default,
                                    "counts": self.Counts,
                                    "txt": self.VariantDetails,
                                    "longtxt": self.LongVariantDetails,
                                    "xls": self.VariantDetailsXLS,
                                    "longxls": self.LongVariantDetailsXLS,
                                    "bed":self.VariantBed,
                                    "featXsamp": self.FeatureXSample,
                                    "featmutXsamp": self.Feature_and_MutationXSample,
                                    "vcf": self.VCF,
                                    "all": self.All,
                                    "runinfo": self.RunInfo }

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
        self.config = config
        self.supported_formats[format]()
        
    def RunInfo(self):
        '''
        Print useful information about the run, including the version, time, and configuration used 

        '''
        outputDirName = self.outputDirName
        config = self.config
        try: 
            ofRunInfo = open(outputDirName + "/run_info.txt", 'w+')
        except:
            abortWithMessage("Error opening output files in {0}/".format(outputDirName))
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
        outputDirName = self.outputDirName
        varDF = self.data
        try:
            ofCounts = open(outputDirName + "/counts.txt", 'w+')
        except:
            abortWithMessage("Error opening output files in {0}/".format(outputDirName))

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

        out.sort().to_csv(ofCounts, sep='\t', na_rep='?', index=True)
        print("\t{0}: {1} rows".format(ofCounts.name, len(out)))        
        return True 

    def FeatureXSample(self):
        '''
        Print counts per feature per sample.
        Sample names populate the table header and feature names populate the first column. The count per sample per feature are the table values.

        Output: feature_by_sample.txt
        '''
        outputDirName = self.outputDirName
        varDF = self.data
        try:
            ofFeatureXSample = pd.ExcelWriter(str(outputDirName) + "/feature_by_sample.xls")
        except:
            abortWithMessage("Error opening output files in {0}/".format(outputDirName))

        groupedDF = pd.DataFrame(varDF.groupby(['feature','sample']).apply(len))
        outDF = groupedDF.stack().unstack(1)
        outDF.index = outDF.index.droplevel(1)
        # check for xls filetype row limit. variant data frames with more than 65,535 lines will casue an error when dumping the data frame.
        if len(outDF) >= 65536:
            throwWarning("featXsamp: There are too many mutations for an Excel xls file. {0} mutations, 65,536 lines maximum.".format(len(outDF)))
            return True
        outDF.to_excel(ofFeatureXSample, 'Feature by Sample', na_rep=0, index=True)
        ofFeatureXSample.save()
        print("\t{0}: {1} rows".format(str(outputDirName) + "/feature_by_sample.xls", len(outDF)))        
        return True 

    def Feature_and_MutationXSample(self):
        '''
        Print counts per mutation per sample.
        Sample names populate the table header and feature names populate the first column. Chromosome, position, ref, and alt populate the next columns.
        The table values are boolean: 1 for present mutation, 0 for missing mutation. 

        Output: feature_and_mutation_by_sample.txt
        '''
        outputDirName = self.outputDirName
        varDF = self.data
        try:
            ofFeature_and_MutationXSample = pd.ExcelWriter(str(outputDirName) + "/feature_and_mutation_by_sample.xls")
        except:
            abortWithMessage("Error opening output files in {0}/".format(outputDirName))

        groupedDF = pd.DataFrame(varDF.groupby(['feature','chr','pos','ref','alt','sample']).apply(len))
        outDF = groupedDF.stack().unstack(5)
        outDF.index = outDF.index.droplevel(5)
        # check for xls filetype row limit. variant data frames with more than 65,535 lines will casue an error when dumping the data frame.
        if len(outDF) >= 65536:
            throwWarning("featmutXsamp: There are too many mutations for an Excel xls file. {0} mutations, 65,536 lines maximum.".format(len(outDF)))
            return True
        outDF.to_excel(ofFeature_and_MutationXSample, 'Feature and Mutation by Sample', na_rep=0, index=True)
        ofFeature_and_MutationXSample.save()
        print("\t{0}: {1} rows".format(str(outputDirName) + "/feature_and_mutation_by_sample.xls", len(outDF)))        
        return True 

    def VariantDetails(self):
        '''
        Print all information about each mutation,
        combining all mutations (irrespective of in how many samples they appear)
        into a single row
        Note: chrom, position, ref, alt, and feature are all required to uniquely identify a mutation 
              indels may have the same chr, pos, but different ref/alt\

        Output: variant_details.txt
        '''
        outputDirName = self.outputDirName
        varDF = self.data

        try:
            ofVariantDetails = open(outputDirName + "/variant_details.txt", 'w+')
        except:
            abortWithMessage("Error opening output files in {0}/".format(outputDirName))
        # Group by (chr, pos, ref, alt, feature)
        grouped = varDF.groupby(['chr', 'pos', 'ref', 'alt', 'feature'])
        # apply collapsing function to each pandas group
        out = grouped.apply(collapseVariantDetails)
        # print the new, collapsed dataframe to a file
        out.sort(['feature','pos']).to_csv(ofVariantDetails, sep='\t', na_rep='?', index=False)
        print("\t{0}: {1} rows".format(ofVariantDetails.name, len(out)))
        return True

    def LongVariantDetails(self):
        '''
        Similar to printVariantDetails above, but writes each instance of a mutation to a new row. 
        Each mutation is written once per source instead of combining reoccurring mutations in to 1 unique row.
        
        Output: long_variant_details.txt
        '''
        outputDirName = self.outputDirName
        varDF = self.data
        try:
            ofLongVariantDetails = open(outputDirName + "/long_variant_details.txt", 'w+')
        except:
            abortWithMessage("Error opening output files in {0}/".format(outputDirName))

        varDF.sort(['feature','pos']).to_csv(ofLongVariantDetails, sep='\t', na_rep='?', index=False)
        print("\t{0}: {1} rows".format(str(outputDirName + '/long_variant_details.txt'), len(varDF)))
        return True

        ofLongVariantDetails.close()
        print("\t{0}: {1} rows".format(ofLongVariantDetails.name, nrow))
        return True

    def VariantDetailsXLS(self):
        '''
        Function to print all information about each mutation,
        combining all mutations (irrespective of in how many samples they appear)
        into a single row
        Note: chrom, position, ref, alt, and feature are all required to uniquely identify a mutation 
              indels may have the same chr, pos, but different ref/alt
        
        Output: variant_details.xls
        Note: switching the pandas ExcelWriter file extension to xlsx instead of xls requires openpyxl
        '''
        outputDirName = self.outputDirName
        varDF = self.data

        try:
            ofVariantDetails = pd.ExcelWriter(str(outputDirName) + '/variant_details.xls')
        except:
            abortWithMessage("Error opening output files in {0}/".format(outputDirName))
        # Group by (chr, pos, ref, alt, feature)
        grouped = varDF.groupby(['chr', 'pos', 'ref', 'alt', 'feature'])
        # apply collapsing function to each pandas group
        out = grouped.apply(collapseVariantDetails)
        # check for xls filetype row limit. variant data frames with more than 65,535 lines will casue an error when dumping the data frame.
        if len(out) >= 65536:
            throwWarning("xls: There are too many mutations for an Excel xls file. {0} mutations, 65,536 lines maximum.".format(len(out)))
            return True
        # print the new, collapsed dataframe to file a
        out.sort(['feature','pos']).to_excel(ofVariantDetails, 'Variant Details', na_rep='?', index=False)
        ofVariantDetails.save()
        print("\t{0}: {1} rows".format(str(outputDirName + '/variant_details.xls'), len(out)))
        return True

    def LongVariantDetailsXLS(self):
        '''
        Similar to printVariantDetails above, but writes each instance of a mutation to a new row. 
        Each mutation is written once per source instead of combining reoccurring mutations in to 1 unique row.
        
        Output: long_variant_details.xls
        Note: switching the pandas ExcelWriter file extension to xlsx instead of xls requires openpyxl
        '''
        outputDirName = self.outputDirName
        varDF = self.data

        # check for xls filetype row limit. variant data frames with more than 65,535 lines will casue an error when dumping the data frame.
        if len(varDF) >= 65536:
            throwWarning("longxls: There are too many mutations for an Excel xls file. {0} mutations, 65,536 lines maximum.".format(len(varDF)))
            return True
        try:
            ofLongVariantDetails = pd.ExcelWriter(str(outputDirName) + '/long_variant_details.xls')
        except:
            abortWithMessage("Error opening output files in {0}/".format(outputDirName))
        varDF.sort(['feature','pos']).to_excel(ofLongVariantDetails, 'Long Variant Details', na_rep='?', index=False)
        ofLongVariantDetails.save()
        print("\t{0}: {1} rows".format(str(outputDirName + '/long_variant_details.xls'), len(varDF)))
        return True

    def VariantBed(self):
        '''
        Print bed file of the variant locations

        Output: variant_locations.bed
        '''
        outputDirName = self.outputDirName
        varDF = self.data
        try:
            ofVariantBeds = open(outputDirName + "/variant_locations.bed", 'w+')
        except:
            abortWithMessage("Error opening output files in {0}/".format(outputDirName))
        grouped = varDF.groupby(['chr', 'pos', 'ref', 'alt'])
        out = grouped.apply(collapseVariantBed)
        out.sort(['chr','start']).to_csv(ofVariantBeds, sep='\t', na_rep='?', index=False, header=False)
        ofVariantBeds.close()
        print("\t{0}: {1} rows".format(ofVariantBeds.name, len(out)))
        return True

    def VCF(self):
        '''
        Print vcf file of the variant locations, feature, depth, and variant frequency.

        Output: variant_locations.vcf
        '''
        outputDirName = self.outputDirName
        varDF = self.data
        try:
            ofVariantVCF = open(outputDirName + "/variant_locations.vcf", 'w+')
        except:
            abortWithMessage("Error opening output files in {0}/".format(outputDirName))
        
        #VCF header to stream
        ofVariantVCF.write("##fileformat=VCFv4.1\n")
        
        grouped = varDF.groupby(['chr', 'pos', 'ref', 'alt'])
        out = grouped.apply(collapseVCF)
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
        self.RunInfo()
        self.VariantDetails()
        self.Counts()

    def All(self):
        '''
        Runs all available output functions, based on those in the self.supported_formats dictionary 
        '''
        formats = self.supported_formats.keys()
        formats.remove('default') # the default output is redundant in this situation, since we are going to run all the functions anwyay.
        formats.remove('all')     # remove 'all' or initiate an infinite loop
        for format in formats:
            self.supported_formats[format]()


############################################
#     V     support functions      V       #
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
            outvals.append( [group[column].unique()[0]] )   # functional consequence and effect should only be 1 value as well, but may not be if there is a mixture of annotated and unannotated vcf files. 
                                                            # This sorts the effects and consequences, and only takes the first value as part of the output. This extracts the annotation and leaves the blank out of the output
                                                            # each of the appended values must be an array in order to transform dictionary to pandas data frame in the next line
    outDF = pd.DataFrame( dict(zip(columns,outvals)), columns=columns )
    return outDF

def collapseVariantBed(group):
    '''
    Pandas operation to support the VariantBed function.
    Collapses variant rows that share the same contig, position, ref allele, and alt allele
    Input: a pandas groupby object
    Output: a pandas dataframe object, ready to be written to bed format
    '''
    bed_fields = ['chr','start','end','name']
    chrom = group['chr'].unique()
    pos = group['pos'].unique()
    start = pos - 1     # bed files are 0-based positions
    end = pos
    name = ";".join(set(group['feature'].values))   # mutations that fall in overlapping features will have a name column with multiple values delimited by semicolons
    return pd.DataFrame( dict( zip(     \
        ['chr','start','end','name'],   \
        [chrom, start , end , name]     \
        )), columns=bed_fields)

def collapseVCF(group):
    '''
    Pandas operation to support the VCF function.
    Collapses variant rows that share the same contig, position, ref allele, and alt allele
    Input: a pandas groupby object
    Output: a pandas dataframe object, ready to be written to vcf format 
    '''
    chrom = group['chr'].unique()
    pos = group['pos'].unique()
    ref = group['ref'].unique()
    alt = group['alt'].unique()
    none = np.array(['.'], dtype=object) # numpy array object for blank columns, such as ID, QUAL, and FILTER
    info = '' # output will have sample-specific information in the 'info' column
    vcf_fields = ['#CHROM','POS','ID', 'REF','ALT', 'QUAL', 'FILTER','INFO']  # standard VCF header
    samples = group['sample'].unique()
    for sample in samples:
        info += str(sample) + ':'
        for feat in group[group['sample'] == sample]['feature'].values:
            line = group[ (group.sample == sample) & (group.feature == feat) ]
            vf = str(line['vf'].values[0])
            dp = str(line['dp'].values[0])
            info += str(feat) + " vf=" + vf + " dp=" + dp + "; "
    return pd.DataFrame(dict(zip(vcf_fields,[chrom, pos, none, ref, alt, none, none, info])), columns=vcf_fields)


def abortWithMessage(message):
    print("*** FATAL ERROR: " + message + " ***")
    exit(2)

def throwWarning(message, help = False):
    print("*** WARNING: " + message + " ***")
    return
