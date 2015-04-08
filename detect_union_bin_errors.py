import argparse
import HTSeq
import pdb # pdb.set_trace()
import itertools
import os
from collections import defaultdict

def FindProblems(featureDict):
    multipleCopies = defaultdict(list)  # keys are gene, values are lists of affected contigs 
                                        # note that not all contigs need to be excluded for an affected gene
    for feature, locs in featureDict.items():
        for contig in locs.keys():
            if len(locs[contig]) > 1:
                # this feature has multiple copies on this contig
                if ThereAreCopies(locs[contig]):
                    # these copies are problematic, ie: they are not overlapping features of the same name
                    multipleCopies[feature].append(contig)
                else:
                    # although there are copies of this feature on the same contig, they overlap each other enough to be 1 bin
                    pass
    return multipleCopies

def ThereAreCopies(locs):
    # if the lowest ending point is less than the highest starting point, these are separate copies
    locLowestEnd = float('inf')
    locHighestStart = 0
    for loc in locs:
        if loc[1] < locLowestEnd:
            locLowestEnd = loc[1]
        if loc[0] > locHighestStart:
            locHighestStart = loc[0]
    if locLowestEnd < locHighestStart:
        return True
    return False

def WriteProblems(problems, outputDir):
    output = open(outputDir + '/union_incompatible_genes.txt', 'w')
    for gene in problems.keys():
        output.write(str(gene + '\n'))
        for contig in problems[gene]:
            output.write(str(gene + '.' + contig + '\n'))
    return True

def main():
    '''
    This tool was designed to aid users running Vaggregate with --union while using feature type = gene_name. 
    It will find features with the same name on the same contig, which cause large, problematic feature bins when running Vaggregate with --union.
    The output is a text document to be placed into the current working directory where Vaggregate will be executed
        It will be read automatically if present at runtimme
    '''
    
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_directory", required=True, help="Directory in which to place output")
    parser.add_argument("-g", "--gff", required=True, help="Annotation GFF/GTF for feature binning")
    parser.add_argument("-f", "--featuretype", required=True, help="Feature type into which to bin. Gencode GTF example: gene_name, gene_id, transcript_name, transcript_id, etc. ")
    args = parser.parse_args()

    gffFile = HTSeq.GFF_Reader(args.gff)
    features = defaultdict(dict)
    for feature in itertools.islice(gffFile, 0, None):
        if feature.type == "gene":
            name = feature.attr[args.featuretype]
            try:
                features[name][feature.iv.chrom].append((feature.iv.start,feature.iv.end))
            except KeyError:
                features[name][feature.iv.chrom] = [(feature.iv.start,feature.iv.end)]

    problems = FindProblems(features)
    WriteProblems(problems, os.path.expanduser(args.output_directory))

if __name__ == "__main__":
    main()