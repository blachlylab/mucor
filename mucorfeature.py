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
#    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

# mucorfeature.py

import HTSeq

# mucor modules
from variant import Variant

class MucorFeature(HTSeq.GenomicFeature):
	'''Specific Genomic Feature. For example, gene SF3B1, band 13q13.1, or chromosome X'''

	def __init__(self, name, type_, interval):
		if name == '': raise NameError('name was an empty string')
		if type_ == '': raise NameError('type_ was an empty string')
		if not isinstance(interval, HTSeq.GenomicInterval): raise TypeError('interval must be of type HTSeq.GenomicInterval')
		self.variants = set()					# empty set to be filled with objects of class Variant
		HTSeq.GenomicFeature.__init__(self, name, type_, interval)
	
	def numVariants(self):
		return len(self.variants)

	def weightedVariants(self):
		'''Instead of returning the number of variants, return the sum of tumor_f for all variants'''
		tumor_f_sum = 0.0
		for var in self.variants:
			tumor_f_sum += float(var.frac)

		return tumor_f_sum

	def uniqueVariants(self):
		'''Return the set of unique variants from the set of all variants (for this feature)'''
		# exploit the hashtable and uniqueness of sets to quickly find
		# unique tuples (contig, pos, ref, alt) of variant info
		# sorted by chrom, pos
		uniqueVariantsTemp = set()
		for var in self.variants:
			candidate = (var.pos.chrom, var.pos.pos, var.ref, var.alt)
			uniqueVariantsTemp.add(candidate)
		# sort by chr, then position
		# TO DO: python sorted() will sort as: chr1, chr10, chr2, chr20, chrX. Fix.
		uniqueVariantsTemp = sorted(uniqueVariantsTemp, key=lambda varx: ( varx[0] + str(varx[1]) ) )

		# Now construct a returnable set of Variant objects,
		# specifying multiple "sources" in the source field
		# this loop's inner-product is #unique variants * #total variants, times #features
		# and is a major inefficiency
		uniqueVariants = set()
		for uniqueVarTup in uniqueVariantsTemp:
			source = ""
			frac = ""   
			dp = ""     
			eff = ""
			fc = ""
			#annot = ""
			for varClass in self.variants:
				if (varClass.pos.chrom, varClass.pos.pos, varClass.ref, varClass.alt) == uniqueVarTup:
					source += varClass.source + ", "
					frac += str(varClass.frac) + ", "   
					dp += str(varClass.dp) + ", "       
					eff += str(varClass.eff) + ", "     
					fc += str(varClass.fc) + ", "
					#annot += str(varClass.annot) + ", " 
			pos = HTSeq.GenomicPosition(uniqueVarTup[0], uniqueVarTup[1] )
			uniqueVar = Variant(source.strip(", "), pos, ref=uniqueVarTup[2], alt=uniqueVarTup[3], frac=str(frac).strip(", "), dp=str(dp).strip(", "), eff=str(eff).strip(", "), fc=str(fc).strip(", ")) ######## Karl Modified ##############
			uniqueVariants.add(uniqueVar)

		return uniqueVariants

	def numUniqueVariants(self):
		'''Return the number of unique variants from the set of all variants (for this feature)'''
		return len(self.uniqueVariants())

	def numUniqueSamples(self):
		sources = set()
		for var in self.variants:
			sources.add(var.source)
		return len(sources)
