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
import sys
try:
    import tabix
except:
    pass

class KnownVariant:
    '''Data about known variants '''
    def __init__(self,source,chrom,pos,ref,alt,rs):
        self.source = source    # database which reports this variant
        self.chrom = chrom      # chromosome/contig
        self.pos = pos          # position
        self.ref = ref          # reference allele 
        self.alt = alt          # alternate allele
        self.rs = rs            # RS number

def abortWithMessage(message):
    print("*** FATAL ERROR: " + message + " ***")
    exit(2)

def isAnnotatedSNP(var, dbs):
    '''
    Check if the input mutation exists in any of the input databases
    Returns a boolean true/false, as well as a string of rs numbers
        If the mutation does not exist in any database, the returned rs number is "?"
        Better than '' when printing output
    '''
    chrom = str(var.pos.chrom)
    '''
    old int tabix: 
    try:
        chrom = int(str(var.pos.chrom).strip('chr'))
    except:
        if str(var.pos.chrom).strip('chr') == "X":
            chrom = int(23)
        elif str(var.pos.chrom).strip('chr') == "Y":
            chrom = int(24)
    '''
    if 'tabix' not in sys.modules:
        return False, str('?')
    try:
        chrom
    except:
        # chromosome is undefined; not a number, not an X or a Y. 
        #     This includes chrM and alternative contigs
        return False, str('?')
        pass
    spos = int(var.pos.pos - 1)     # start position
    epos = int(var.pos.pos)         # end position
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
                for row in tb.query(chrom, spos, epos):
                    if str(row[3]) == ref and str(row[4]) in str(alt).split(','):
                        datab += str(row[2]) + ", "
    if str(datab) != "":
        return True, str(datab.strip(', '))
    elif str(datab) == "":
        return False, str('?')
