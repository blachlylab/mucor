class Info:
    '''Program info: logo, version, and usage'''
    logo = """
 __   __   ______     ______     ______     ______     ______     ______     ______     ______   ______    
/\ \ / /  /\  __ \   /\  ___\   /\  ___\   /\  == \   /\  ___\   /\  ___\   /\  __ \   /\__  _\ /\  ___\   
\ \ \\'/   \ \  __ \  \ \ \__ \  \ \ \__ \  \ \  __<   \ \  __\   \ \ \__ \  \ \  __ \  \/_/\ \/ \ \  __\   
 \ \__|    \ \_\ \_\  \ \_____\  \ \_____\  \ \_\ \_\  \ \_____\  \ \_____\  \ \_\ \_\    \ \_\  \ \_____\ 
  \/_/      \/_/\/_/   \/_____/   \/_____/   \/_/ /_/   \/_____/   \/_____/   \/_/\/_/     \/_/   \/_____/ 
                                                                                                                                          
"""
    version = "0.9"
    versionInfo = "vaggregate version {0}\nJames S Blachly, MD\nKarl W Kroll, BS".format(version)
    # usage needs to be updated, or eliminated if argparse can replace this function. 
    usage = """
Usage:
{0} [-h] | -g featurefile.gff -f feature_type [-u] -o <outputDir> <mutect001.txt mutect002.txt ... mutectNNN.txt>

Flags:
    -h  Show this help

    -g  GFF3/GTF file describing the features into which mutations will be binned

    -f  String to describe the type of feature to bin for this run.
        e.g. gene_id or transcript_id or chromosome_id

    -u,
    --union
        Join all items with same ID for feature_type (specified by -f)
        into a single, continuous bin. For example, if you want intronic
        variants counted with a gene, use this option. 
        ** TO DO **
        WARNING, this will likely lead to spurious results due to things
        like MIR4283-2 which exists twice on + and - strand of the same
        chromosome over 1 megabase apart. This creates one huge spurious
        bin.

    -o  Specify output directory

Output directory:
    Output files will be placed in the specified directory.
    If the directory already exists, an error will occur (won't overwrite)

    Output consists of CSV spreadsheets (.txt):
        1. Plain text report of # variants binned according to feature_type
        2. Summary of pre-filter and post-filter variants
        3. Detailed report of all variants by feature_type

Input files:
    <mutect001.txt mutect002.txt ... mutectNNN.txt>
    Final arguments should be a list of mutations in muTect output format

"""
    description = """

vaggregate: Variant Aggregation

vaggregate reads in variant files from a variety of sources (VCF, muTect .out)
and counts the number of mutations falling into known features. These are
grouped together and output to see which features (genes) and which spec-
-ific locations within those genes have the highest frequency of mutation
within a group.
"""
    epilog = """

Output directory:
    Output files will be placed in the specified directory.
    If the directory already exists, an error will occur (won't overwrite)

    Output consists of CSV spreadsheets (.txt):
        1. Plain text report of # variants binned according to feature_type
        2. Summary of pre-filter and post-filter variants
        3. Detailed report of all variants by feature_type

Input files:
    <mutect001.txt mutect002.txt ... mutectNNN.txt>
    Final arguments should be a list of mutations in muTect output format
"""
