# databases.py
#
# define genomic databases such as:
# dbSNP
# 1000 Genomes
# 6k Exomes
# COSMIC


class databases
# makes a dictionary out of dbSNP,
# with tuple of chrom,position as the key; and the rs number as values.
#
def load_dbsnp():
	startTime = time.clock()
	snps = defaultdict(str)
	dbsnp_p = '/nfs/17/osu7366/projects/new_AK/dbSNPandMiSeq.P'
	#dbsnp_file = '/nfs/17/osu7366/reference/snp138Common.txt.gz'
	#dbsnp = gzip.open(dbsnp_file,'rb')
	print("\n=== Reading dbSNP pickle file {0} ===".format(dbsnp_p))
	'''
	for line in dbsnp:
		col = line.split('\t')
		if str(col[11]) == "deletion": # deletions in our VCF file start 1 base upstream (-1) from dbSNP, but have the correct rs number
			snps[tuple((str(col[1]), int(col[3]) - 1))] = str(col[4])
		else:
			snps[tuple((str(col[1]), int(col[3])))] = str(col[4])
	'''
	snps = pickle.load(open(dbsnp_p,'rb'))
	totalTime = time.clock() - startTime
	print("{0} sec\t{1} SNPs".format(int(totalTime), len(snps.values())))
	return snps

##################################

######## Karl Added ##############
# true or false to check if a location
# (tuple of chrom,position) is in the dbSNP dictionary. 
# must use defaultdict above to avoid key errors here
def in_dbsnp(snps, loc):
	status = False
	annotation = snps[loc]
	if str(annotation).startswith('rs'):
		status = True
	return status
##################################
