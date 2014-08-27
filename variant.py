# variant.py

class Variant:
	'''Data about SNV and Indels'''
	def __init__(self,source,pos,ref,alt,frac,dp, eff):
		self.source = source	# source of variant - typically a filename
		self.pos = pos			# HTSeq.GenomicPosition
		self.ref = ref
		self.alt = alt
		self.frac = frac        ######## Karl Added ##############
		self.dp = dp            ######## Karl Added ##############
		self.eff = eff          ######## Karl Added ##############
		#self.annot = annot      ######## Karl Added ##############
