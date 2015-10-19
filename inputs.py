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

# inputs.py
from __future__ import print_function
import HTSeq
import sys
from variant import Variant
from pdb import set_trace as stop

def throwWarning(message, help = False):
    print("*** WARNING: " + message + " ***")
    return

class Parser(object):
    '''Object to cover all parsing functions'''

    def __init__(self):
            self.row = ''
            self.source = ''
            self.var = None
            self.fieldId = ''
            self.header = ''
            self.fn  = ''

            # Respective parse functions. 
            # Keys are identical to the allowed formats in mucor_config.py
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

    def parse(self, row, source):
        self.source = source
        self.row  = row
        chrom = row.chrom 
        ref   = row.ref
        alt   = "/".join(sorted(row.alt))
        pos   = row.pos # already exists as GenomicPosition object
        dp    = 0
        frac  = 0.0
        effect, fc = parse_EFC(row.info)
        out = []
        samples_dict = self.supported_formats[self.source](row.samples)
        for sample, vals in samples_dict.items():
            out.append(Variant(source='', sample=sample, pos=pos, ref=ref, alt=alt, frac=vals[1], dp=vals[0], eff=effect, fc=fc))
        return out

    def parse_MiSeq(self,samples):
        ''' MiSeq vcf parser function. Input: InputParser object. Output: Variant object '''
        out = {} # list of variant objects; 1 per sample
        for sample, values in samples.items():
            vf = None
            dp = None 
            if 'VF' in values.keys():
                vf = float(values['VF'])
            if 'DP' in values.keys():
                dp = int(values['DP'])

            if 'AD' in values.keys():
                # prepare ref depth and alt depth in case VAF or DP was not defined 
                rd = int(values['AD'].split(',')[0])
                ad = int(values['AD'].split(',')[1])
            if not vf and 'AD' in values.keys():
                # try to calculate vaf from AD column
                vf = float(ad)/float(ad + rd)
            if not dp and 'AD' in values.keys():
                # try to calculate dp from AD column
                if vf:
                    # if AD and VF are both defined, this will yield the same DP defined in the INFO column
                    dp = round(ad/vf)
                else:
                    # the sum of ref and alt alleles may not equal total depth, if there are two mutations at the same location, or if a single base was mutated and not shown
                    dp = rd + ad
            out[sample] = (dp, vf) 
        return out

    def parse_IonTorrent(self, samples):
        ''' Ion Torrent vcf parser function. Input: InputParser object. Output: Variant object '''
        out = {} # list of variant objects; 1 per sample
        for sample, values in samples.items():
            dp = None 
            vf = None 
            try:
                ao = sum( [int(x) for x in values['AO'].split(',')] )
                ro = int(values['RO'])
                dp = int(values['DP'])
                vf = float(ao)/float(dp)
            except KeyError:
                pass
            out[sample] = (dp, vf) 
        return out

    def parse_MuTectOUT(self):
        ''' MuTect '.out' parser function. Input: InputParser object. Output: Variant object '''

        row = self.row
        fieldId = self.fieldId
        header = self.header
        fn = self.fn 
        chrom = row[0]
        ref   = row[3]
        alt   = row[4]
        effect = ""
        fc = ""
        vf = float(row[fieldId['tumor_f']])
        dp = int(int(str(row[fieldId['t_ref_count']]).strip()) + int(str(row[fieldId['t_alt_count']]).strip()))
        position = int(row[fieldId['position']])

        var = Variant(source=fn.split('/')[-1], sample='', pos=HTSeq.GenomicPosition(chrom, int(position)), ref=ref, alt=alt, frac=vf, dp=dp, eff=effect.strip(';'), fc=fc.strip(';'))
        return var
        
    def parse_MuTectVCF(self, samples):
        ''' MuTect vcf parser function. Input: InputParser object. Output: Variant object '''
        out = {} # list of variant objects; 1 per sample
        for sample, values in samples.items():
            dp = None
            vf = None 
            # remove empty sample columns
            if sample == "none" and values['DP'] == '0' and values['BQ'] == '.':
                continue
            try:
                dp = values['DP']
                vf = values['FA']
            except KeyError:
                pass
            out[sample] = (dp, vf) 
        return out

    def parse_SomaticIndelDetector(self, samples):
        ''' DEPRECIATED AND UNTESTED '''
        ''' GATK SomaticIndelDetector vcf parser function. Input: InputParser object. Output: Variant object '''
        out = {} # list of variant objects; 1 per sample
        for sample, values in samples.items():
            dp = None
            vf = None 
            try:
                dp = int(values['DP'])
                ad = int(values['AD'].split(',')[1])
                vf = float(ad)/float(dp)
            except KeyError:
                pass
            out[sample] = (dp, vf) 
        return out

        '''
        # pre-htseq vcf reder parse function
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
        '''

    def parse_SamTools(self, samples):
        ''' samtools vcf parser function. Input: InputParser object. Output: Variant object '''
        out = {} # list of variant objects; 1 per sample
        for sample, values in samples.items():
            dp = None
            vf = None 
            if 'DP4' in values.keys():
                dp4 = values['DP4'].split(',')
                ro = int(dp4[0]) + int(dp4[1])
                ao = int(dp4[2]) + int(dp4[3])
                dp = ro + ao
                vf = float(ao)/float(dp)
            if not dp:
                # if DP4 is not present, check for depth field
                if 'DP' in values.keys():
                    dp = int(values['DP'])

            out[sample] = (dp, vf) 
        return out

    def parse_VarScan(self,samples):
        ''' varscan vcf parser function. Input: InputParser object. Output: Variant object '''
        out = {} # list of variant objects; 1 per sample
        for sample, values in samples.items():
            dp = None
            vf = None 
            try:
                dp = int(values['DP'])
                vf = float(values['FREQ'].strip('%'))/100.0
            except KeyError:
                pass
            out[sample] = (dp, vf) 
        return out

    def parse_HapCaller(self, samples):
        ''' GATK haplotype caller vcf parser function. Input: InputParser object. Output: Variant object '''
        out = {} # list of variant objects; 1 per sample
        for sample, values in samples.items():
            vf = None
            dp = None 
            try:
                dp = int(values['DP'])
                ad = str(values['AD'])
                ref_count = int(ad.split(',')[0])
                alt_count = int(ad.split(',')[1])
                vf = float( float(alt_count)/(float(ref_count) + float(alt_count)) )
            except KeyError:
                pass
            out[sample] = (dp, vf)
        return out

    def parse_FreeBayes(self, samples):
        ''' freebayes vcf parser function. Input: InputParser object. Output: Variant object '''
        out = {} # list of variant objects; 1 per sample
        for sample, values in samples.items():
            dp = None
            vf = None
            try:
                dp = int(values['DP'])
                ro = int(values['RO'])
                ao = sum([int(x) for x in values['AO'].split(',')])
                vf = float( float(ao)/float(ao + ro) )
            except  KeyError:
                pass
            out[sample] = (dp, vf) 
        return out

    def parse_GenericGATK(self, samples):
        ''' 
        Generic GATK parser function. This was written for the Illumina BaseSpace BWA Enrichment Workflow vcf files, but may apply to more filetypes
        Input: InputParser object. Output: Variant object 
        '''
        out = {}
        for sample, values in samples.items():
            dp = None
            vf = None
            if 'DP' in values.keys():
                #DP is defined 
                dp = values['DP']
            if 'VF' in values.keys():
                #VF is defined
                vf = values['VF']
            if not dp or not vf:
                try:
                    # try to determine depth and VAF from 'AD' field, if present
                    ro = int(values['AD'].split(',')[0])
                    ao = sum([int(x) for x in values['AD'].split(',')[1:]])
                    if not dp:
                        dp = ro + ao
                    vf = float(float(ao)/float(dp)) # one VF for all possible alternate alleles. Nothing unusual, unless the mutation has multiple alt alleles in 1 vcf line
                except KeyError:
                    if not dp:
                        throwWarning("Could not determine depth for {0} using GenericGATK parser".format(sample))
                    if not vf:
                        throwWarning("Could not determine variant allele frequency for {0} using GenericGATK parser".format(sample))
                    pass
            out[sample] = (dp, vf) 
        return out

    def parse_MAF(self):
        ''' DEPRECIATED AND UNTESTED'''
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

def parse_EFC(INFO):
    # attempt to extract 'effect' and 'functional consequence' from the VCF line
    effect = ""
    fc = ""
    muts = []
    loca = []
    try:
        for eff in INFO.split(';'):
            if eff.startswith('EFF='):
                for j in eff.split(','):
                    muts.append(str(j.split('|')[3]))
                    loca.append(str(j.split('(')[0]).replace('EFF=',''))
        # reformat the lists to exclude blanks and be semicolon delimited
        for mut in set(muts):
            if str(mut) != "":
                effect += str(mut) + ";"
        for loc in set(loca):
            if str(loc) != "":
                fc += str(loc) + ";"
    except KeyError:
        # this VCF may not be snpEff annotated
        pass

    # if the MiSeq software reported functional consequence and effect and the file is not snpEff anotated, the MiSeq annotations will be used instead
    for i in INFO.split(';'):
        if i.startswith("FC=") and not fc:
            for j in i.split('=')[1].split(','):
                if str(j.split('_')[0]) not in str(fc):
                    fc += str(j.split('_')[0]) + ";"
                try:
                    if str(j.split('_')[1]) not in str(effect):
                        effect += str(j.split('_')[1]) + ";"
                except IndexError:
                    pass
        elif str(i) == "EXON":
            fc += 'EXON'

    return effect.strip(';'), fc.strip(';')