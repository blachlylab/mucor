# report.py
#
# create, format, and write an executive report
# on the muCorrelated sample variants

import os
import pandas as pd

class Report(object):
    '''Create, format, and write an executive report on the muCorrelated sample variants'''
    
    reporttypes = ['latex']    
    templates = ['report.tex']
    
    def __init__(self, study, variantsDataFrame):
        if study == '': raise ValueError('study name was an empty string')
        if not isinstance(variantsDataFrame, pd.DataFrame): raise TypeError('variantsDataFrame was not an instance of pandas.DataFrame')
        self.study = study
        self.vardf = variantsDataFrame

    def write(self, filename, filetype="latex"):
        if filename == '': raise ValueError('filename was an empty string')
        if filetype not in self.reporttypes: return ValueError('filetype {} unrecognized'.format(filetype))
        pass

    @classmethod
    def check_templates(cls):
        '''Ensure the report templates are available. Raise exception if not.'''
        dirname = os.path.dirname(os.path.abspath(__file__))
        for templ in cls.templates:
            if not os.path.isfile(os.path.join(dirname,templ)):
                raise RuntimeError('Report template {} not found'.format(os.path.join(dirname,templ)))

