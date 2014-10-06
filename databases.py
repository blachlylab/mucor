# databases.py
#
# define genomic databases such as:
# dbSNP
# 1000 Genomes
# 6k Exomes
# COSMIC
import gzip
import pdb
from collections import defaultdict
import time

class KnownVariant:
	'''Data about known variants '''
	def __init__(self,source,chrom,pos,ref,alt,rs):
		self.source = source	# database which reports this variant
		self.chrom = chrom		# chromosome/contig
		self.pos = pos			# position
		self.ref = ref 			# reference allele 
		self.alt = alt			# alternate allele
		self.rs = rs			# RS number

def abortWithMessage(message):
    print("*** FATAL ERROR: " + message + " ***")
    exit(2)

def str_to_bool(s):
    if str(s) == 'False':
        return False
    else:
        return True

def parseDB(db, snps, duplicateAnnots):  # parse database in VCF format 
	if str(db).split('.')[-1] == str("gz"):
		database = gzip.open(db)
	elif str(db).split('.')[-1] == str("vcf"):
		database = open(db)
	else: abortWithMessage("Error opening database files: {0}".format(db))
	try:
		row = database.readline()
	except StopIteration: 
		print("Empty file {}".format(db))
	while str(row)[0:2] == '##':
		if str("source=") in str(row):
			source = str(row).split("=")[1].strip()
		row = database.readline()
	header = row
	if len(header) == 0: 
		raise ValueError('Invalid header')
	fieldId = dict(zip(header, range(0, len(header))))
	for line in database:
		col = line.split('\t')
		chrom = col[0]
		pos   = col[1]
		rs 	  = col[2]
		ref   = col[3]
		alt   = col[4]
		#source= str("dbSNP")
		if not snps.has_key((chrom,pos)):
			snps[(chrom,pos)] = KnownVariant(source=[source],chrom=chrom,pos=pos,ref=ref,alt=alt,rs=[rs])
		elif snps.has_key((chrom,pos)): 
			if snps[(chrom,pos)].ref == ref and snps[(chrom,pos)].alt == alt:
				if snps[(chrom,pos)].rs != rs:
					snps[(chrom,pos)].rs.append(rs)
					snps[(chrom,pos)].source.append(source)
				elif snps[(chrom,pos)].rs == rs:
					pass
			else:
				duplicateAnnots.add(snps[(chrom,pos)])
				duplicateAnnots.add(KnownVariant(source=[source],chrom=chrom,pos=pos,ref=ref,alt=alt,rs=[rs]))
				'''
				pdb.set_trace()
				abortWithMessage("Error: {0} from {1} already exists as a different SNV".format( str((chrom,pos)), str(db).split('/')[-1] ))
				'''

	return snps


def load_db(dbs):
	startTime = time.clock()
	snps = {}
	duplicateAnnots = set()
	for db in dbs:
		if not str_to_bool(db):
			pass
		else:
			print str(db)
			print("\n=== Reading db file {0} ===".format(db))
			parseDB(db, snps, duplicateAnnots)
			totalTime = time.clock() - startTime
			print("{0} sec\t{1} SNPs".format(int(totalTime), len(snps.values())))
	if duplicateAnnots:
		print("*** WARNING: {0} mutations have identical positions, but different REF/ALTs".format(len(duplicateAnnots)))
	return snps

##################################

######## Karl Added ##############
# true or false to check if a location
# (tuple of chrom,position) is in the dbSNP dictionary. 
# must use defaultdict above to avoid key errors here
def isAnnotatedSNP(snps, loc):
	if snps.has_key(loc):
		return bool(True)
	elif not snps.has_key(loc):
		return bool(False)

##################################
