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
    dbVAFs = {}

    for source, db in dbs.items():
        if not os.path.exists(db):
            pass
        else:
            # 'source' is a variable used to title the column in the output
            # it is defined by the user in the configuration script step when generating the JSON file
            # TO-DO: make a database object to perform and document these file checks. Goal is to run them once per database, rather than once per mutation. 
            #        will improve runtime
            if os.path.splitext(db)[1] == ".gz" and os.path.exists(db + ".tbi"):
                try:
                    database = gzip.open(db)
                except:
                    print("WARNING: could not open {}".format(db))
                    continue
            elif os.path.splitext(db)[1] == ".vcf":
                abortWithMessage("Error: database file {0} must compressed with bgzip".format(db))
            elif os.path.splitext(db)[1] == ".gz" and not os.path.exists(db + ".tbi"):
                abortWithMessage("Compressed database is not tabix indexed")
            else: abortWithMessage("Error opening database files: {0}".format(db))
            
            try:
                row = database.readline()
            except StopIteration: 
                print("Empty file {}".format(db))

            # perform the database query 
            tb = tabix.open(db)
            if len(str(var.alt).split(',')) >= 1:
                dbEntries[source] = '?'
                dbVAFs[source] = '?'
                try:
                    for row in tb.query(var.pos.chrom, spos, epos):
                        if str(row[3]) == var.ref and str(row[4]) in str(var.alt).split(','):
                            # 3rd column (zero indexed = [2])
                            # in a VCF is the ID
                            ID = row[2]

                            # check for known population allele frequencies in the database VCF
                            if "AF=" in row[7]:
                                for item in row[7].split(';'):
                                    if item.startswith('AF='):
                                        dbVAFs[source] = item.split('=')[1]
                                        break

                            if ID == ".":
                                # a '.' indicates it is annotated, but no RS number. 
                                # different from a '?', which indicates not-annotated
                                # this is more clear
                                ID = "PRESENT"
                            dbEntries[source] = ID
                            
                except tabix.TabixError:
                    pass
                    #there are no mutations on this variant's contig, in this particular database 
                        
    return dbEntries, dbVAFs
