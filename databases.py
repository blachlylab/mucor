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
import tabix

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
			snps[(chrom,pos)] = ([source],str(chrom),int(pos),str(ref),str(alt),[rs]) #{'source':[source],'chrom':str(chrom),'pos':int(pos),'ref':str(ref),'alt':str(alt),'rs':[rs]}
		elif snps.has_key((chrom,pos)): 
			if snps[(chrom,pos)][3] == ref and snps[(chrom,pos)][4] == alt:
				if snps[(chrom,pos)][5] != rs:
					snps[(chrom,pos)][5].append(rs)
					snps[(chrom,pos)][0].append(source)
				elif snps[(chrom,pos)][5] == rs:
					pass
			else:
				duplicateAnnots.append(snps[(chrom,pos)])
				duplicateAnnots.append(([source],str(chrom),int(pos),str(ref),str(alt),[rs]))
				'''
				pdb.set_trace()
				abortWithMessage("Error: {0} from {1} already exists as a different SNV".format( str((chrom,pos)), str(db).split('/')[-1] ))
				'''

	return snps


def load_db(dbs):
	startTime = time.clock()
	snps = {}
	duplicateAnnots = list()
	for db in dbs:
		if not bool(db):
			pass
		else:
			print("\n=== Reading db file {0} ===".format(db))
			parseDB(db, snps, duplicateAnnots)
	totalTime = time.clock() - startTime
	if bool(db):
		print("{0} sec\t{1} SNPs".format(int(totalTime), len(snps.values())))
	if duplicateAnnots:
		print("*** WARNING: {0} mutations have identical positions, but different REF/ALTs".format(len(duplicateAnnots)))
	return snps

######## Karl Added ##############
# true or false to check if a location
# (tuple of chrom,position) is in the dbSNP dictionary. 
# must use defaultdict above to avoid key errors here
'''
def isAnnotatedSNP(snps, loc):
	if snps.has_key(tuple(loc)):
		return bool(True)
	elif not snps.has_key(loc):
		return bool(False)
'''
def isAnnotatedSNP(var, dbs):
	try:
		chrom = int(str(var.pos.chrom).strip('chr'))
	except:
		if str(var.pos.chrom).strip('chr') == "X":
			chrom = int(23)
		elif str(var.pos.chrom).strip('chr') == "Y":
			chrom = int(24)
	try:
		chrom
	except:
		# chromosome is undefined; not a number, not an X or a Y. 
		#     This includes chrM and alternative contigs
		return False, str('?')
		pass
	spos = int(var.pos.pos - 1)
	epos = int(var.pos.pos)
	ref = var.ref
	alt = var.alt
	datab = ""
	for db in dbs:
		if not bool(db):
			pass
		else:
			if str(db).split('.')[-1] == str("gz"):
				database = gzip.open(db)
			elif str(db).split('.')[-1] == str("vcf"):
				abortWithMessage("Error: database file {0} must compressed with bgzip".format(db))
			else: abortWithMessage("Error opening database files: {0}".format(db))
			try:
				row = database.readline()
			except StopIteration: 
				print("Empty file {}".format(db))
			while str(row)[0:2] == '##':
				if str("source=") in str(row):
					source = str(row).split("=")[1].strip()
				row = database.readline()
			tb = tabix.open(db)
			if len(str(alt).split(',')) >= 1:
				for row in tb.queryi(chrom, spos, epos):
					if str(row[3]) == ref and str(row[4]) in str(alt).split(','):
						datab += str(row[2]) + ", "
	if str(datab) != "":
		return True, str(datab.strip(', '))
	elif str(datab) == "":
		return False, str('?')

##################################
