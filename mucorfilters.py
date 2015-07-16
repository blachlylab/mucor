# mucorfilters.py
#

class MucorFilters(object):
	"""docstring for MucorFilters"""
	def __init__(self):
		super(MucorFilters, self).__init__()
		self._filterSets = None					# internally use dict
		self._filters = None					# internally use dict

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

