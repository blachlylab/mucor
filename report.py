# report.py
#
# create, format, and write an executive report
# on the muCorrelated sample variants

import pandas as pd

class Report(object):
    '''Create, format, and write an executive report on the muCorrelated sample variants'''

    def __init__(self, study, variantsDataFrame):
        if study == '': raise ValueError('study name was an empty string')
        if not isinstance(variantsDataFrame, pd.DataFrame): raise TypeError('variantsDataFrame was not an instance of pandas.DataFrame')
        self.study = study
        self.vardf = variantsDataFrame

    def write(self, filename, filetype="latex"):
        if filename == '': raise ValueError('filename was an empty string')
        if filetype not in ['latex']: return ValueError('filetype {} unrecognized'.format(filetype))
        pass

