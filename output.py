#
# output.py
#
# James S Blachly, MD
# 2014-10-04

import os
import pandas as pd

class Writer(object):
    """Object that parses the mucor dataframe and can write output in several different formats"""
    
    def __init__(self):
        self.data = pd.DataFrame()
        
        #self.supported_formats = ["default", "long"] # FUTURE: VCF, GVF, XLS
        self.supported_formats = {  "default": self.default,
                                    "long": self.long }

    def write(self, data, format, outdir):
        """Write data in format to outdir
        data: a pandas dataframe
        format: a string: one of [default, vcf, gvf, xls]   # consider merging in report functionality
        outdir: a string"""
        
        # Parameter integrity/validity checks
        if (type(data) == pd.DataFrame):
            self.data = data   
        else:
            raise TypeError("data was not pandas.core.frame.DataFrame")
        
        if (format not in self.supported_formats):
            raise ValueError("Unsupported output format {}".format(format)) 
        
        if not os.path.exists(outdir):
            raise ValueError("The specified output directory ({}) does not exist.".format(outdir))

        supported_formats[format]()


    def default():
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

