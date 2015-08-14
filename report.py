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

# report.py
#
# create, format, and write an executive report
# on the aggregated sample variants

import os
import pandas as pd

class Report(object):
    '''Create, format, and write an executive report on the aggregated sample variants'''
    
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

