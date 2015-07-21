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

# variant.py

class Variant:
	'''Data about SNV and Indels'''
	def __init__(self,source,pos,ref,alt,frac,dp, eff,fc):
		self.source = source	# source of variant - typically a filename
		self.pos = pos			# HTSeq.GenomicPosition
		self.ref = ref
		self.alt = alt
		self.frac = frac     
		self.dp = dp            
		self.eff = eff          
		self.fc = fc

