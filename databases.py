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

# databases.py
#
# define genomic databases such as:
# dbSNP
# 1000 Genomes
# 6k Exomes
# COSMIC
import gzip
from collections import defaultdict
import time
import os
import sys
try:
    import tabix
except ImportError:
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

def checkAndOpen(db):
    if 'tabix' not in sys.modules:
        return None

    db = os.path.expanduser(db)
    if not os.path.exists(db):
        pass
    else:
        # 'source' is a variable used to title the column in the output
        # it is defined by the user in the configuration script step when generating the JSON file
        if os.path.splitext(db)[1] == ".gz" and os.path.exists(db + ".tbi"):
            try:
                database = gzip.open(db)
            except IOError:
                print("WARNING: could not open {}".format(db))
                return None
        elif os.path.splitext(db)[1] == ".vcf":
            abortWithMessage("Error: database file {0} must compressed with bgzip".format(db))
        elif os.path.splitext(db)[1] == ".gz" and not os.path.exists(db + ".tbi"):
            abortWithMessage("Compressed database is not tabix indexed")
        else: abortWithMessage("Error opening database files: {0}".format(db))
        
        try:
            row = database.readline()
        except StopIteration: 
            print("Empty file {}".format(db))
            return None

        tb = tabix.open(db)
        return tb

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

    spos = int(var.pos.pos - 1)     # start position
    epos = int(var.pos.pos)         # end position
    alts = str(var.alt).split(',')  # list of 1 or more alternative alleles 
    dbEntries = {}

    for source, tb in dbs.items():
        if len(alts) >= 1:
            dbEntries[source] = '?'
            dbEntries[source + "_VAF"] = '?'
            try:
                for row in tb.query(var.pos.chrom, spos, epos):
                    if str(row[3]) == var.ref and str(row[4]) in alts:
                        # 3rd column (zero indexed = [2])
                        # in a VCF is the ID
                        ID = row[2]

                        # check for known population allele frequencies in the database VCF
                        if "AF=" in row[7]:
                            for item in row[7].split(';'):
                                if item.startswith('AF='):
                                    dbEntries[source + "_VAF"] = float(item.split('=')[1])
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
                        
    return dbEntries
