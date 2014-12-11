# inputs.py
import HTSeq
from variant import Variant

class InputParser(object):
    '''Object to cover all parsing functions'''

    def __init__(self):
            self.row = ''
            self.fieldId = ''
            self.header = ''
            self.fn = ''
            self.eff = ''
            self.fc = ''

def parse_MiSeq(InputParser):
    ''' MiSeq vcf parser function. Input: InputParser object. Output: Variant object '''

    row = InputParser.row
    fieldId = InputParser.fieldId
    header = InputParser.header
    fn = InputParser.fn 
    chrom = row[0]
    ref   = row[3]
    alt   = row[4]
    fc    = InputParser.fc
    effect   = InputParser.eff
    for i in row[fieldId['INFO']].split(';'):
        if i.startswith("DP="):
            dp = i.split('=')[1]
        if i.startswith("FC=") and not fc:
            global SnpEff_switch
            SnpEff_switch = bool(True)
            for j in i.split('=')[1].split(','):
                if str(j.split('_')[0]) not in str(fc):
                    fc += str(j.split('_')[0]) + ";"
                try:
                    if str(j.split('_')[1]) not in str(effect):
                        effect += str(j.split('_')[1]) + ";"
                except:
                    pass
        elif str(i) == "EXON":
            fc += 'EXON'
    if not fc:
        fc = str("?")
    if not effect:
        effect = str("?")
    k = 0
    for i in row[fieldId['FORMAT']].split(':'):
        if str(i) == "VF":
            vf = float( row[fieldId[header[-1]]].split(':')[k] )
        '''
        #for when vf is not in the format column, but AD is
        if str(i) == "AD" and not dp or not vf:
            dp = 0
            rd = int(row[fieldId[header[-1]]].split(':')[k].split(',')[0])
            ad = int(row[fieldId[header[-1]]].split(':')[k].split(',')[1])
            dp = int(rd) + int(ad)
        '''
        k += 1

    position = int(row[fieldId['POS']])
    #return vf, dp, position, eff, fc
    var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
    return var

def parse_IonTorrent(InputParser):
    ''' Ion Torrent vcf parser function. Input: InputParser object. Output: Variant object '''

    row = InputParser.row    
    fieldId = InputParser.fieldId
    header = InputParser.header
    fn = InputParser.fn 
    chrom = row[0]
    ref   = row[3]
    alt   = row[4]
    effect = InputParser.eff
    fc = InputParser.fc
    for i in row[fieldId['INFO']].split(';'):
        if i.startswith("AO="):
            tempval = i.split('=')[1]
        if i.startswith("RO="):
            ro = i.split('=')[1]
        if i.startswith("DP="):
            dp = i.split("=")[1]
    if str(',') in str(tempval):
        tempval2 = [int(numeric_string) for numeric_string in tempval.split(',')]
        try:
            ao = sum(tempval2)
        except:
            abortWithMessage("AO should be an int, or a list of ints: AO = {0}/".format(tempval2))
    else:
        ao = tempval
    vf = float(float(ao)/float(float(ro) + float(ao)))
    position = int(row[fieldId['POS']])
    for i in str(row[fieldId['ALT']]).split(','):
        if len(str(row[fieldId['REF']])) > len(i):
            #this is a deletion in Ion Torrent data
            position = int(row[fieldId['POS']])
            break
    #return vf, dp, position
    var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
    return var

def parse_MuTectOUT(InputParser):
    ''' MuTect '.out' parser function. Input: InputParser object. Output: Variant object '''

    row = InputParser.row
    fieldId = InputParser.fieldId
    header = InputParser.header
    fn = InputParser.fn 
    chrom = row[0]
    ref   = row[3]
    alt   = row[4]
    effect = InputParser.eff
    fc = InputParser.fc
    vf = row[fieldId['tumor_f']]
    dp = int(int(str(row[fieldId['t_ref_count']]).strip()) + int(str(row[fieldId['t_alt_count']]).strip()))
    position = int(row[fieldId['position']])

    #return vf, dp, position 
    var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
    return var
    
def parse_MuTectVCF(InputParser):
    ''' MuTect vcf parser function. Input: InputParser object. Output: Variant object '''

    row = InputParser.row
    fieldId = InputParser.fieldId
    header = InputParser.header
    fn = InputParser.fn 
    chrom = row[0]
    ref   = row[3]
    alt   = row[4]
    effect = InputParser.eff
    fc = InputParser.fc
    j = 0
    for i in InputParser.header:
        if str(i) not in ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'none']: # This line should detect if the sample id is in the line.
            tmpsampID = i
    for i in row[fieldId['FORMAT']].split(':'):
        if i == "FA":
            vf = row[fieldId[tmpsampID]].split(':')[j]
        elif i == "DP":
            dp = row[fieldId[tmpsampID]].split(':')[j]
        j+=1
    position = int(row[fieldId['POS']])

    #return vf, dp, position 
    var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
    return var

def parse_SomaticIndelDetector(InputParser):
    ''' GATK SomaticIndelDetector vcf parser function. Input: InputParser object. Output: Variant object '''

    row = InputParser.row
    fieldId = InputParser.fieldId
    header = InputParser.header
    fn = InputParser.fn 
    chrom = row[0]
    ref   = row[3]
    alt   = row[4]
    effect = InputParser.eff
    fc = InputParser.fc
    j = 0
    # Below attempts to grab sample ID. 
    # assumes that sample ID is the final column in the InputParser.header. always true? 
    # if not always true, adopt the parse_mutect solution here as well
    tmpsampID = header[-1]
    
    for i in row[fieldId['FORMAT']].split(':'):
        if i == "AD":
            ALT_count = row[fieldId[tmpsampID]].split(':')[j].split(',')[1]
        elif i == "DP":
            dp = row[fieldId[tmpsampID]].split(':')[j]
            vf = float( float(ALT_count)/float(dp) )
        j+=1
    position = int(row[fieldId['POS']])
    #return vf, dp, position
    var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
    return var

def parse_SamTools(InputParser):
    ''' samtools vcf parser function. Input: InputParser object. Output: Variant object '''

    row = InputParser.row
    fieldId = InputParser.fieldId
    header = InputParser.header
    fn = InputParser.fn 
    chrom = row[0]
    ref   = row[3]
    alt   = row[4]
    effect = InputParser.eff
    fc = InputParser.fc
    position = int(row[fieldId['POS']])
    for i in row[fieldId['INFO']].split(';'):
        if i.startswith("DP4="):
            j = i.split('=')[1].split(',')
            ref = int(int(j[0]) + int(j[1]))
            alt = int(int(j[2]) + int(j[3]))
            dp = int(int(ref) + int(alt))
            vf = float( float(alt)/float(dp) )
            #return vf, dp, position
            var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
            return var

def parse_VarScan(InputParser):
    ''' varscan vcf parser function. Input: InputParser object. Output: Variant object '''

    row = InputParser.row
    fieldId = InputParser.fieldId
    header = InputParser.header
    fn = InputParser.fn 
    chrom = row[0]
    ref   = row[3]
    alt   = row[4]
    effect = InputParser.eff
    fc = InputParser.fc
    j = 0
    position = int(row[fieldId['POS']])
    for i in row[fieldId['FORMAT']].split(':'):
        if str(i) == "DP":
            dp = int(row[fieldId[header[-1]]].split(':')[j])
        if str(i) == "FREQ":
            vf = float(float(str(row[fieldId[header[-1]]].split(':')[j]).strip('%'))/float(100))
        j += 1
    #return vf, dp, position
    var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
    return var

def parse_HapCaller(InputParser):
    ''' GATK haplotype caller vcf parser function. Input: InputParser object. Output: Variant object '''

    row = InputParser.row
    fieldId = InputParser.fieldId
    header = InputParser.header
    fn = InputParser.fn 
    chrom = row[0]
    ref   = row[3]
    alt   = row[4]
    effect = InputParser.eff
    fc = InputParser.fc
    j = 0
    position = int(row[fieldId['POS']])
    for i in row[fieldId['FORMAT']].split(':'):
        if str(i) == "DP":
            dp = int(row[fieldId[header[-1]]].split(':')[j])
        if str(i) == "AD":
            ad = str(row[fieldId[header[-1]]].split(':')[j])
            if str(',') in ad:
                ref = int(ad.split(',')[0])
                alt = int(ad.split(',')[1])
                vf = float( float(alt)/(float(ref) + float(alt)) )
            else:
                abortWithMessage("Sample {0} may not have Haplotype Caller mutations with no ALT or vf".format(header[-1]))
        j += 1
    #return vf, dp, position
    var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
    return var

def parse_FreeBayes(InputParser):
    ''' freebayes vcf parser function. Input: InputParser object. Output: Variant object '''

    row = InputParser.row
    fieldId = InputParser.fieldId
    header = InputParser.header
    fn = InputParser.fn 
    chrom = row[0]
    ref   = row[3]
    alt   = row[4]
    effect = InputParser.eff
    fc = InputParser.fc
    j = 0
    position = int(row[fieldId['POS']])
    for i in row[fieldId['FORMAT']].split(':'):
        if str(i) == "DP":
            dp = int(row[fieldId[header[-1]]].split(':')[j])
        if str(i) == "RO":
            ro = int( str(row[fieldId[header[-1]]].split(':')[j]) )
        if str(i) == "AO":
            ao = int(  sum([ int(x) for x in str(row[fieldId[header[-1]]].split(':')[j]).split(',')]) )
        j += 1
    vf = float( float(ao)/float(ao + ro) )
    #return vf, dp, position
    var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
    return var

def parse_MAF(InputParser):
    ''' maf filetype parser function. Input: InputParser object. Output: Variant object '''

    row = InputParser.row
    fieldId = InputParser.fieldId
    header = InputParser.header
    fn = InputParser.fn 
    position = int(str(row[fieldId['Start_position']]).split('.')[0]) # case sensitive. what if, 'Start_Position' instead? case-insensitive hash lookup, or make everything lowercase befor making comparisons?
    dp = int(str(row[fieldId['TTotCov']]).split('.')[0])
    vf = float( float(row[fieldId['TVarCov']])/float(dp) )
    chrom = str(row[fieldId['Chromosome']])
    ref =  str(row[fieldId['Reference_Allele']])
    alt = str(row[fieldId['Tumor_Seq_Allele2']])
    effect = InputParser.eff
    fc = InputParser.fc
    if ref == "-":
        ref = ""
    if alt == "-":
        alt = ""
    #return vf, dp, position, chrom, ref, alt    
    var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
    return var
