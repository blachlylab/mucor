# inputs.py
import HTSeq
from variant import Variant
import pdb

class Parser(object):
    '''Object to cover all parsing functions'''

    def __init__(self):
            self.row = ''
            self.fieldId = ''
            self.header = ''
            self.fn = ''
            self.eff = ''
            self.fc = ''
            self.source = ''

            # Respective parse functions. 
            # Keys are identical to the allowed formats in vagg_config.py
            # Values are the parse functions defined below
            self.supported_formats = {  "MiSeq":self.parse_MiSeq,
                                        "IonTorrent":self.parse_IonTorrent,
                                        "SomaticIndelDetector":self.parse_SomaticIndelDetector,
                                        "Mutect":self.parse_MuTectVCF,
                                        "muTector":self.parse_MuTectOUT,
                                        "Samtools":self.parse_SamTools,
                                        "VarScan":self.parse_VarScan,
                                        "HaplotypeCaller":self.parse_HapCaller,
                                        "FreeBayes":self.parse_FreeBayes,
                                        "GenericGATK":self.parse_GenericGATK }

    def parse(self, source, row, fieldId, header, fn, effect, fc):
        self.source = source
        self.row  = row
        self.fieldId = fieldId
        self.header = header
        self.fn  = fn
        self.eff  = effect
        self.fc = fc
        return self.supported_formats[self.source]()

    def parse_MiSeq(self):
        ''' MiSeq vcf parser function. Input: InputParser object. Output: Variant object '''

        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        fc    = self.fc
        effect   = self.eff
        for i in row[fieldId['INFO']].split(';'):
            if i.startswith("DP="):
                dp = i.split('=')[1]

            # if the MiSeq software reported functional consequence and effect and the file is not snpEff anotated, the MiSeq annotations will be used instead
            if i.startswith("FC=") and not fc:
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
        var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var

    def parse_IonTorrent(self):
        ''' Ion Torrent vcf parser function. Input: InputParser object. Output: Variant object '''

        row = self.row    
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        effect = self.eff
        fc = self.fc
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
                # this is a deletion in Ion Torrent data
                position = int(row[fieldId['POS']])
                break
        var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var

    def parse_MuTectOUT(self):
        ''' MuTect '.out' parser function. Input: InputParser object. Output: Variant object '''

        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        effect = self.eff
        fc = self.fc
        vf = float(row[fieldId['tumor_f']])
        dp = int(int(str(row[fieldId['t_ref_count']]).strip()) + int(str(row[fieldId['t_alt_count']]).strip()))
        position = int(row[fieldId['position']])

        var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var
        
    def parse_MuTectVCF(self):
        ''' MuTect vcf parser function. Input: InputParser object. Output: Variant object '''

        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        effect = self.eff
        fc = self.fc
        j = 0
        for i in self.header:
            if str(i) not in ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'none']: # This line should detect if the sample id is in the line.
                tmpsampID = i
        for i in row[fieldId['FORMAT']].split(':'):
            if i == "FA":
                vf = float(row[fieldId[tmpsampID]].split(':')[j])
            elif i == "DP":
                dp = row[fieldId[tmpsampID]].split(':')[j]
            j+=1
        position = int(row[fieldId['POS']])

        var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var

    def parse_SomaticIndelDetector(self):
        ''' GATK SomaticIndelDetector vcf parser function. Input: InputParser object. Output: Variant object '''

        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        effect = self.eff
        fc = self.fc
        j = 0
        # Below attempts to grab sample ID. 
        # assumes that sample ID is the final column in the self.header. always true? 
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
        var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var

    def parse_SamTools(self):
        ''' samtools vcf parser function. Input: InputParser object. Output: Variant object '''

        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        effect = self.eff
        fc = self.fc
        position = int(row[fieldId['POS']])
        for i in row[fieldId['INFO']].split(';'):
            if i.startswith("DP4="):
                j = i.split('=')[1].split(',')
                ro = int(int(j[0]) + int(j[1]))
                ao = int(int(j[2]) + int(j[3]))
                dp = int(int(ro) + int(ao))
                vf = float( float(ao)/float(dp) )
                var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
                return var

    def parse_VarScan(self):
        ''' varscan vcf parser function. Input: InputParser object. Output: Variant object '''

        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        effect = self.eff
        fc = self.fc
        j = 0
        position = int(row[fieldId['POS']])
        for i in row[fieldId['FORMAT']].split(':'):
            if str(i) == "DP":
                dp = int(row[fieldId[header[-1]]].split(':')[j])
            if str(i) == "FREQ":
                vf = float(float(str(row[fieldId[header[-1]]].split(':')[j]).strip('%'))/float(100))
            j += 1
        var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var

    def parse_HapCaller(self):
        ''' GATK haplotype caller vcf parser function. Input: InputParser object. Output: Variant object '''

        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        effect = self.eff
        fc = self.fc
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
        var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var

    def parse_FreeBayes(self):
        ''' freebayes vcf parser function. Input: InputParser object. Output: Variant object '''

        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        effect = self.eff
        fc = self.fc
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
        var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var

    def parse_GenericGATK(self):
        ''' 
        Generic GATK parser function. This was written for the Illumina BaseSpace BWA Enrichment Workflow vcf files, but may apply to more filetypes
        Input: InputParser object. Output: Variant object 
        '''
        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        effect = self.eff
        fc = self.fc
        j = 0
        position = int(row[fieldId['POS']])
        for i in row[fieldId['FORMAT']].split(':'):
            if str(i) == "AD":
                ro = int(row[fieldId[header[-1]]].split(':')[j].split(',')[0])
                ao = int(row[fieldId[header[-1]]].split(':')[j].split(',')[1])
                dp = ro + ao
                vf = float(float(ao)/float(dp))
                break
            j += 1

        var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var

    def parse_MAF(self):
        ''' maf filetype parser function. Input: InputParser object. Output: Variant object '''

        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        position = int(str(row[fieldId['Start_position']]).split('.')[0]) # case sensitive. what if, 'Start_Position' instead? case-insensitive hash lookup, or make everything lowercase befor making comparisons?
        dp = int(str(row[fieldId['TTotCov']]).split('.')[0])
        vf = float( float(row[fieldId['TVarCov']])/float(dp) )
        chrom = str(row[fieldId['Chromosome']])
        ref =  str(row[fieldId['Reference_Allele']])
        alt = str(row[fieldId['Tumor_Seq_Allele2']])
        effect = self.eff
        fc = self.fc
        if ref == "-":
            ref = ""
        if alt == "-":
            alt = ""   
        var = Variant(source=fn.split('/')[-1], pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var
