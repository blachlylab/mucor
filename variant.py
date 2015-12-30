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

# variant.py

from __future__ import print_function

class Variant:
	'''Data about SNV and Indels'''
	def __init__(self,source,sample,pos,ref,alt,frac,dp, eff,fc):
		self.source = source	# source of variant - typically a filename
		self.sample = sample	# sample ID, as defined in the VCF 
		self.pos = pos			# HTSeq.GenomicPosition
		self.ref = ref
		self.alt = alt
		self.frac = frac     
		self.dp = dp            
		self.eff = eff          
		self.fc = fc
	def __str__(self):
		out = ""
		for k,v in {'source':self.source, 'sample':self.sample, 'pos':self.pos, 'ref':self.ref, 'alt':self.alt, 'frac':self.frac, 'dp':self.dp, 'eff':self.eff, 'fc':self.fc}.items():
			out += k + ":\t" + str(v) + "\n"
		return out.strip()
