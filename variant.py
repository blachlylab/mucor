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

