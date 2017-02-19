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

# mucorfilters.py
#
from collections import defaultdict

class MucorFilters(object):
    """docstring for MucorFilters"""
    def __init__(self):
        super(MucorFilters, self).__init__()
        self._filterSets = None					# internally use dict
        self._filters = None					# internally use dict
        self.vcfFilters = []
        self.regionDict = defaultdict(set)

	def loadFromXML(self, rootNode):
		'''Load filters and filterSets from mucor XML config file.
		Overwrites any existing filters or filterSets stored in this object.'''

		# namespace
		ns = "{urn:osuccc:hcrg:mucor}"

		if not isinstance(rootNode, ET.Element):
			raise TypeError
		if root.tag != ns + "mucor":
			raise ValueError

		fNode = root.find(ns + "filters")
		fsNode = root.find(ns + "filterSets")
		if fNode == None or fsNode == None:
			raise ValueError

		self._filters = dict()
		self._filterSets = dict()

		for mcFilter in fNode:		# 'filter' is a reserved keyword
			if mcFilter.tag != ns + 'filter':
				raise ValueError
			if 'id' not in mcFilter.attrib.keys():
				raise ValueError

			# add empty string entry to the _filters dictionary
			filterName = mcFilter.attrib['id']
			self._filters[filterName] = ""

			# build the filter expression
			fieldId = mcFilter.find(ns + 'fieldId')
			comparison = mcFilter.find(ns + 'comparison')
			value = mcFilter.find(ns + 'value')

			if fieldId == None or comparison == None or value == None:
				raise ValueError
			if fieldId.text == None or comparison.text == None or value.text == None:
				raise ValueError

			comparison_dict = {
				'equal' : '==',
				'notequal' : '!=',
				'lessthan' : '<',
				'lessorequal' : '<=',
				'greaterthan' : '>',
				'greaterorequal' : '>='
			}

			filter_expr = '( row.' + fieldId.text + ' ' \
						+ comparison_dict[comparison.text] + " '" \
						+ value.text + "' )"

			self._filters[filterName] = filter_expr

		for filterSet in fsNode:
			if filterSet.tag != ns + 'filterSet':
				raise ValueError
			if 'id' not in filterSet.attrib.keys():
				raise ValueError

			# add empty list entry to the _filterSets dictionary
			filterSetName = filterSet.attrib['id']
			self._filterSets[filterSetName] = list()

			# populate the list
			for filterId in filterSet.findall(ns + 'filterId'):	# if no filterId tags found, this is ok, no error (an empty filterSet)
				if filterId.text == None:						# however if a filterId is not specified, this is malformed XML
					raise ValueError

				filterName = filterId.text
				filterExpr = self.filter(filterName)
				assert filterExpr != None

				self._filterSets[filterSetName].append(filterExpr)


	def filter(self, filterName):
		'''Return the python expr (as string) corresponding to filterName, or None if filterName not in the _filters dictionary'''
		if self._filters == None: return None
		elif filterName not in self._filters.keys():
			return None
		else:
			return self._filters[filterName]	# return the python expr as string

	@property
	def filters(self):
		'''Return a list of filterIds (names)'''
		if self._filters == None: return None
		else:
			return self._filters.keys()

	@property
	def filterSets(self):
		'''Return a list of filterSet ids (names)'''
		if self._filterSets == None: return None
		else:
			return self._filterSets.keys()

	def filtersInSet(self, filterSetName):
		'''Return a list of python evaluable statements representing the filters in the given filterSet'''
		if self._filterSets == None: return None
		if filterSetName not in self._filterSets.keys(): return None

		return self._filterSets[filterSetName]
    
    def parseRegionBed(regionfile, regionDictIn):
        ''' 
        Read through the supplied bed file and add the rows to the input region dictionary before returning it. Appends to the input dict, rather than overwriting it.
        '''
        regionDict = regionDictIn
        for line in open(str(regionfile),'r'):
            col = line.strip().split("\t")
            chrom = col[0]
            start = col[1]
            end = col[2]
            if len(col) > 3:
                name = col[3]
            else:
                name = str(chrom) + ":" + str(start) + "-" + str(end)
            regionDict[chrom].add((start, end, name))
        return regionDict

    def generateRegionFilter(self, regions):
        for item in regions:
            if os.path.splitext(item)[1].lower() == 'bed':      # this item is a bed file
                self.regionDict = parseRegionBed(item, regionDict)
            elif str(str(item).split(':')[0]).startswith('chr'):    # this is a string 
                reg_chr = str(item.split(':')[0])
                try:
                    reg_str = str(str(item.split(':')[1]).split('-')[0])
                    reg_end = str(str(item.split(':')[1]).split('-')[1])
                except IndexError:                                  # represent whole chromosome regions [ex: chr2] by chrN:0-0 in the region dictionary   
                    reg_str = 0
                    reg_end = 0
                regionDict[reg_chr].add((reg_str, reg_end))
        return

    def filterLoc(self, chrom, start, end):
        '''
        Checks the given mutation location to see if it is in the dictionary of regions
        returning True means the loc should be filterd out/ignored since it is not in the region dict
        returning False means the loc should not be excluded; it is in the regions of interest.  
        '''
        if not self.regionDict[chrom]: # are there any regions of interest in the same chromosome as this mutation?
            return False
        else:
            for locs in self.regionDict[chrom]: # breaking the list of regions according to chromosome should significantly decrease the number of comparisons necessary 
                if locs[0] == 0 and locs[1] == 0: # chrN:0-0 is used to define an entire chromosome as a region of interest. 
                    return False
                elif int(start) >= int(locs[0]) and int(end) <= int(locs[1]):
                    return False

    def filterVCFRow(self, row, kind, fieldId=None):
        '''
        Checks the vcf filter and/or mutect judgement column for permitted filters
        If a variant data row has multiple filters, they must all be permitted for the row to pass
  
        returning True means this row will be filtered out [masked]
        returning False means the row will not be filtered out [pass filter]
        '''
        if kind in ["vcf", "vcf.gz"]:
            for rowFilter in row.filter.split(';'):    ## VCF file format
                if rowFilter not in self.vcfFilters:
                    return True
            if self.filterLoc(row.pos.chrom, int(row.pos.pos), int(row.pos.pos)):
                return True
        if kind == "out":
            for rowFilter in row[fieldId['judgement']].split(';'): ## MuTect '.out' file format
                if rowFilter not in self.vcfFilters:
                    return True
        return False
