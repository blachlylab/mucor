#
# output.py
#
# James S Blachly, MD
# 2014-10-04

import os
import pandas as pd
from config import Config
import pdb

class Writer(object):
    """Object that parses the mucor dataframe and can write output in several different formats"""
    
    def __init__(self):
        self.data = pd.DataFrame()
        self.config = Config()
        self.outputDirName = ''
        self.supported_formats = {  "runinfo": self.RunInfo,
                                    "default": self.Default,
                                    "counts": self.Counts,
                                    "txt": self.VariantDetails,
                                    "longtxt": self.LongVariantDetails,
                                    "xls": self.VariantDetailsXLS,
                                    "longxls": self.LongVariantDetailsXLS }

    def write(self, data, format, outputDirName, config):
        """Write data in format to outputDirName
        data: a pandas dataframe
        format: a string: one of [default, vcf, gvf, xls]   # consider merging in report functionality
        outputDirName: a string"""
        
        # Parameter integrity/validity checks
        if (type(data) == pd.DataFrame):
            self.data = data   
        else:
            raise TypeError("data was not pandas.core.frame.DataFrame")
        
        if (format not in self.supported_formats):
            raise ValueError("Unsupported output format {}".format(format)) 
        
        if not os.path.exists(outputDirName):
            raise ValueError("The specified output directory ({}) does not exist.".format(outputDirName))
        else:
            self.outputDirName = outputDirName
        self.config = config
        try:
            self.supported_formats[format]()
        except:
            pdb.set_trace()

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
        Replacement function

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
        numHits.name = 'numHits'

        uniqueHits = grouped['sample'].count().groupby(level=0).count()
        uniqueHits.name = 'uniqueHits'

        numSamples = grouped['sample'].nunique()
        numSamples.name = 'numSamples'

        # TO DO: Construct DF from the 3 columns
        # TO DO: Write the DF to ofCounts
        
        return True    

    def VariantDetails(self):
        '''
        Replacement function to print all information about each mutation,
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
        out.sort(['chr','pos','feature']).to_csv(ofVariantDetails, sep='\t', na_rep='?', index=False)
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

        varDF.sort(['chr','pos','feature']).to_csv(ofLongVariantDetails, sep='\t', na_rep='?', index=False)
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
        '''
        outputDirName = self.outputDirName
        varDF = self.data

        try:
            ofVariantDetails = pd.ExcelWriter(str(outputDirName) + '/variant_details.xls')
        except:
            pdb.set_trace()
            #abortWithMessage("Error opening output files in {0}/".format(outputDirName))
        # Group by (chr, pos, ref, alt, feature)
        grouped = varDF.groupby(['chr', 'pos', 'ref', 'alt', 'feature'])
        # apply collapsing function to each pandas group
        out = grouped.apply(collapseVariantDetails)
        # print the new, collapsed dataframe to a file
        out.sort(['chr','pos','feature']).to_excel(ofVariantDetails, 'Variant Details', na_rep='?', index=False)
        ofVariantDetails.save()
        print("\t{0}: {1} rows".format(str(outputDirName + '/variant_details.xls'), len(out)))
        return True

    def LongVariantDetailsXLS(self):
        '''
        Similar to printVariantDetails above, but writes each instance of a mutation to a new row. 
        Each mutation is written once per source instead of combining reoccurring mutations in to 1 unique row.
        
        Output: long_variant_details.xls
        '''
        outputDirName = self.outputDirName
        varDF = self.data

        ofLongVariantDetails = pd.ExcelWriter(str(outputDirName) + '/long_variant_details.xls')
        varDF.sort(['chr','pos','feature']).to_excel(ofLongVariantDetails, 'Long Variant Details', na_rep='?', index=False)
        ofLongVariantDetails.save()
        print("\t{0}: {1} rows".format(str(outputDirName + '/long_variant_details.xls'), len(varDF)))
        return True

##########

    def Default(self):
        print "hello world"
        pdb.set_trace()
        raise ValueError("Output format default not implemented ... sorry")

    def long():
        """Print output format in long (record-based) format, where each row represents a single variant in a single case; this is functionally equivalent to catting and sorting all input files"""
            ## ## ##
            # print record format (long format) all variants data frame
            # TO DO : break out into function to slicing and dicing
            #if (varDF['chr'][0].lower().startswith(('chr', 'Chr'))):
            #    varDF['chr_num'] = varDF['chr'].apply(lambda x:x[3:]).astype(int)   # strip off the 'chr', unstringify
            #    varDF.sort(columns=['chr_num', 'pos'], inplace=True)
            #   varDF.drop(['chr_num'], axis=1, inplace=True)
            #else: 
            
            # sorting by chr does not sort properly: lexical order is chr1,chr10,...,chr2,...
            # so we will sort by feature (gene etc.) then position
        self.data.sort(columns=['feature','pos'], inplace=True)
        self.data.replace('', np.nan, inplace=True)
        self.data.to_csv(self.outputDir + '/allvars.txt', sep='\t', na_rep='?', index=False)

    def vcf():
        """Print output in multi-sample VCF 4.1 format"""
        raise ValueError("Output format vcf is not implemented ... sorry")

    def gvf():
        """Print output in multi-sample GVF format"""
        raise ValueError("Output format gvf is not implemented ... perhaps oneday")

def collapseVariantDetails(group):
    '''
    Pandas operations to support the printVariantDetails family of functions. 
    Collapses variant rows that share the same contig, position, ref allele, alt allele, and feature.
    Input: a pandas groupby object
    Output: a pandas dataframe object
    '''
    outvals = []
    columns = list(group.keys().values)
    for column in columns:
        outstring = ''
        if column in ['vf', 'dp', 'sample', 'source']:  # the only columns that need to be concatenated 
                                                        # the others are uniquified and "always" yield 1 value
            for i in group[column].values:
                outstring += str(i) + ", "
            outvals.append( outstring[:-2] )            # trim the extra ', ' off the end of outstring 
        else:
            outvals.append( group[column].unique() )
    outDF = pd.DataFrame( dict(zip(columns,outvals)), columns=columns )
    return outDF

def abortWithMessage(message):
    print("*** FATAL ERROR: " + message + " ***")
    exit(2)