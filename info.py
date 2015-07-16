class Info:
    '''Program info: logo, version, and usage'''
    logo = """

 __    __     __  __     ______     ______     ______    
/\ "-./  \   /\ \/\ \   /\  ___\   /\  __ \   /\  == \   
\ \ \-./\ \  \ \ \_\ \  \ \ \____  \ \ \/\ \  \ \  __<   
 \ \_\ \ \_\  \ \_____\  \ \_____\  \ \_____\  \ \_\ \_\ 
  \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_/ /_/ 
                                                         

"""
    version = "1.0"
    versionInfo = "mucor version {0}\nJames S Blachly, MD\nKarl W Kroll, BS".format(version)
    # usage needs to be updated, or eliminated if argparse can replace this function. 
    usage = """
Usage:
{0} [-h] [-ex] -g GFF -f FEATURETYPE
                       [-db <dbName:/path/database.vcf.gz>] -s
                       <sample_list.txt> [-d <dirname>] [-vcff VCF_FILTERS]
                       [-a ARCHIVE_DIRECTORY] [-r REGIONS] [-u] -jco
                       JSON_CONFIG_OUTPUT -outd OUTPUT_DIRECTORY
                       [-outt OUTPUT_TYPE]

Flags:
  -h, --help            show this help message and exit
  -ex, --example        Print a valid, example JSON config file and exit.
  -g GFF, --gff GFF     Annotation GFF/GTF for feature binning
  -f FEATURETYPE, --featuretype FEATURETYPE
                        Feature type into which to bin. Gencode GTF example:
                        gene_name, gene_id, transcript_name, transcript_id,
                        etc.
  -db <dbName:/path/database.vcf.gz>, --database <dbName:/path/database.vcf.gz>
                        Colon delimited name and path to variant database in
                        bgzipped VCF format. Can be declared >= 0 times.
  -s <sample_list.txt>, --samples <sample_list.txt>
                        Text file containing sample names. One sample per
                        line.
  -d <dirname>, --project_directory <dirname>
                        Working/project directory, in which to find input
                        variant call files. Default: current working directory
  -vcff VCF_FILTERS, --vcf_filters VCF_FILTERS
                        Comma separated list of VCF filters to allow. Default:
                        PASS
  -a ARCHIVE_DIRECTORY, --archive_directory ARCHIVE_DIRECTORY
                        Specify directory in which to read/write archived
                        annotations. Undeclared will prevent using the
                        annotation archive features.
  -r REGIONS, --regions REGIONS
                        Comma separated list of bed regions and/or bed files
                        by which to limit output. Ex:
                        chr1:10230-10240,chr2,my_regions.bed
  -u, --union           Join all items with same ID for feature_type
                        (specified by -f) into a single, continuous bin. For
                        example, if you want intronic variants counted with a
                        gene, use this option. WARNING, this will lead to
                        spurious results for features that are duplicated on
                        the same contig. When feature names are identical, the
                        bin will range from the beginning of the first
                        instance to the end of the last, even if they are
                        several megabases apart. Refer to the documentation
                        for a resolution using 'detect_union_bin_errors.py'
  -jco JSON_CONFIG_OUTPUT, --json_config_output JSON_CONFIG_OUTPUT
                        Name of JSON configuration output file
  -outd OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Name of directory in which to write Mucor output
  -outt OUTPUT_TYPE, --output_type OUTPUT_TYPE
                        Comma separated list of desired output types. Options
                        include: [longtxt, all, featXsamp, longxls, runinfo,
                        counts, txt, vcf, default, bed, mutXsampVAF, mutXsamp,
                        xls]. Default: counts,txt

"""
    description = """

mucor: Mutation Correlation

mucor reads in variant files from a variety of sources (VCF, muTect .out)
and counts the number of mutations falling into known features. These are
grouped together and output to see which features (genes) and which spec-
-ific locations within those genes have the highest frequency of mutation
within a group.
"""
    epilog = """

Output directory:
    Output files will be placed in the specified directory.
    If the directory already exists, an error will occur (won't overwrite)

"""
