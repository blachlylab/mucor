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
import os
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

def dbLookup(var, dbs):
    '''
    Check if the input mutation exists in any of the input databases
    Returns:
        a dictionary {'database source' => 'entry ID (e.g. rs#)'}

        If the mutation does not exist in any database, the returned id number is "?"
        Better than '' when printing output // TO DO: consider '.' or removal entirely
    '''

    if 'tabix' not in sys.modules:
        return False, str('disabled')

    spos = int(var.pos.pos - 1)     # start position
    epos = int(var.pos.pos)         # end position

    dbEntries = {}

    for db in dbs:
        if not os.path.exists(db):
            pass
        else:
            # 'source' is a variable used to title the column in the output
            # it defaults to the database VCF filename, but if a source= pragma
            # is present in the VCF header, it will use that instead
            source = os.path.split(db)[1]
            if os.path.splitext(db)[1] == ".gz":
                try:
                    database = gzip.open(db)
                except:
                    print("WARNING: could not open {}".format(db))
                    continue
            elif os.path.splitext(db)[1] == ".vcf":
                abortWithMessage("Error: database file {0} must compressed with bgzip".format(db))
            else: abortWithMessage("Error opening database files: {0}".format(db))
            
            try:
                row = database.readline()
            except StopIteration: 
                print("Empty file {}".format(db))

            # read the pragma lines;
            # if the optional source= pragma is present,
            # store it instead of filename
            while str(row)[0:2] == '##':
                if str("source=") in str(row):
                    source = str(row).split("=")[1].strip()
                row = database.readline()

            tb = tabix.open(db)
            if len(str(var.alt).split(',')) >= 1:
                for row in tb.query(chrom, spos, epos):
                    if str(row[3]) == var.ref and str(row[4]) in str(var.alt).split(','):
                        # 3rd column (zero indexed = [2])
                        # in a VCF is the ID
                        dbEntries[source] = row[2]
                    else:
                        dbEntries[source] = '?'     # TO DO: Consider change to '.' or removal entirely

    return dbEntries
