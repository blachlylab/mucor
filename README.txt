Introduction
-----------------------------
Mucor is a software to aggregate mutation information from multiple VCF files into a variety of summary files with varying detail. The outputs range in utility, including high-level summaries and full, detailed reports, useful to both biologists and computation members.

System Requirements
-----------------------------
Python Modules:
standard python modules include os, sys, time, argparse, csv, itertools, collections, gzip, cPickle, and json. Nonstandard modules include numpy, pandas, and HTSeq. Optional modules include tabix.

Annotation Files:
Mucor requires a reference annotation file, such as the GTF/GFF associated with the reference genome of the relevant organism. 
Known mutations databases are optional if the user elects to decorate output with dbSNP, 1000 genome, or other mutation databases. These files must be in bgzipped VCF format and indexed with tabix. Users may construct their own annotation, provided it conforms to the established VCF  standards. 

Supported Variant Formats:
Mucor can parse VCF files created by the Ion Torrent, MiSeq, Somatic Indel Detector, Mutect, Samtools, VarScan, Haplotype Caller, FreeBayes, and GATK softwares. 

Installation
-----------------------------
clone github - TBA

Files and Operation
-----------------------------
User Interface Scripts:
mucor_config.py - Configure the run by creating a JSON file to pass to mucor.py
mucor.py - The main mucor script that parses data and generates output

Mucor Functions and Modules:
config.py - Contains the Config class, used in mucor.py to assign the JSON file to an object
databases.py - Contains a class and function used in mucor.py o check for known, annotated variants
inputs.py - Contains the Parser class, used in mucor.py to read variant call files
mucorfeature.py - Contains the MucorFeature class, used in mucor.py to store feature data
output.py - Contains the Writer class, used in mucor.py to write various output files
variant.py - Contains the Variant class, used in mucor.py to store data about each variant

Operation
-----------------------------
Operating mucor occurs in two distinct steps with their own python scripts; configuration and execution. 

Configuration: mucor_config.py
The configuration step is completed using the provided `mucor_config.py` utility. It accepts the following command line arguments and creates a JSON file that will be passed to the main Mucor script. All required parameters must be supplied. Optional parameters may be left blank, some of which will assume the given default values. 

  -ex, --example        Print a valid, example JSON config file and exit.
Function that will write a template of the JSON config. It can be edited manually and supplied to mucor. 

  -g GFF, --gff GFF
Reference annotation GFF/GTF for feature binning. Required

  -f FEATURETYPE, --featuretype FEATURETYPE 
Feature type into which to bin. Gencode GTF example: gene_name, gene_id, transcript_name, transcript_id, etc. Required

  -db DATABASES, --databases DATABASES
Colon delimited name and path to variant database in bgzipped VCF format. Can be declared >= 0 times. Ex: -db name1:/full/user/path/name1.vcf.gz. Optional

  -s SAMPLES, --samples SAMPLES 
Text file containing sample names. One sample per line. Configuration may not be correct if any sample names are within another sample name. Ex: U-23 and U-238. U-23 would erroneously identify U-238 files, requiring manual modification of the JSON file. Required

  -d PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY 
Working/project directory, in which to find variant call files to aggregate. Variant calls can be in the provided directory, or any of its subdirectories. Default: current working directory

  -vcff VCF_FILTERS, --vcf_filters VCF_FILTERS 
Comma separated list of VCF filters to allow. Default: PASS

 -a ARCHIVE_DIRECTORY, --archive_directory ARCHIVE_DIRECTORY
Specify directory in which to read/write archived annotations. This step will significantly speed up future runs that use the same annotation and feature type, even if the sample data changes. Undefined will prevent using the annotation archive features. Optional

  -r REGIONS, --regions REGIONS 
Comma separated list of bed regions and/or bed files by which to limit output. Bed regions can be specific positions, or entire chromosomes. Ex: chr1:10230-10240,chr2,my_regions.bed. Optional

  -jco JSON_CONFIG_OUTPUT, --json_config_output JSON_CONFIG_OUTPUT
Name of JSON configuration output file. This is the configuration file fed into Mucor. Required

  -outd OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
Name of directory in which to write Mucor output. Required

  -outt OUTPUT_TYPE, --output_type OUTPUT_TYPE
Comma separated list of desired output types. Options include: counts, txt, longtxt, xls, longxls, bed, featXsamp, featmutXsamp, all. Default: counts,txt

Execution: mucor.py
Mucor is executed by running the main, mucor.py python script with the desired json configuration file as the only argument. This facilitates rerunning analyses in the future and documenting the selected parameters to maintain consistency between runs. 

python mucor.py configuration.json

Outputs
-----------------------------
Users may select any number of output formats from the list below. These must be specified in the JSON configuration file prior to executing the main Mucor script.  

counts
Print counts of mutations per feature. Output: counts.txt

txt
Print all information about each mutation, combining all instances of a mutation (irrespective of in how many samples it appears) into a single row. Useful for mutation-centric output for a project. Identical to xls in layout, but different in file format. Output: variant_details.txt

longtxt
Similar to txt above, but writes each instance of a mutation to a new row. Each mutation is written once per source file, instead of combining reoccurring mutations in to 1 unique row. Identical to longxls in layout, but different in file format. Output: long_variant_details.txt

xls
Print all information about each mutation, combining all instances of a mutation (irrespective of in how many samples it appears) into a single row. Useful for mutation-centric output for a project. Identical to txt in layout, but different in file format. Output: variant_details.xls

longxls
Similar to xls above, but writes each instance of a mutation to a new row. Each mutation is written once per source instead of combining reoccurring mutations in to 1 unique row. Identical to longtxt in layout, but different in file format. Output: long_variant_details.xls

bed
Print bed file of the variant locations Output: variant_locations.bed

vcf
Print vcf file of the variant locations, features, depths, and variant frequencies. Output: variant_locations.vcf

featXsamp
Print table of mutation counts per feature per sample. Sample names populate the table header and feature names populate the first column. The count of unique mutations per sample per feature are the table values. Output: feature_by_sample.txt

featmutXsamp
Print table of mutations per sample. Sample names populate the table header and feature names populate the first column. Chromosome, position, ref, and alt populate the second, third, and fourth columns respectively. The table values are boolean: 1 for present mutation, 0 for missing mutation. Output: feature_and_mutation_by_sample.txt

all
Execute all output types

Known Issues
-----------------------------
There is an issue when a variant file presents a contig that the pickled (archived) annotation does not have. The solution is to disable the archive feature by omitting the -a, --archive_directory option. This will permit the unknown contig in output, but all mutations on the unknown contig will be shown as having no feature. It may also be possible to make a new pickle file using an annotation GTF/GFF which does contain the contig in question.

  File "mucor/mucor.py", in parseVariantFiles
    resultSet = gas[ var.pos ]      # returns a set of zero to n IDs (e.g. gene symbols)
  File "_HTSeq.pyx", line 521, in HTSeq._HTSeq.GenomicArray.__getitem__ (src/_HTSeq.c:10700)
KeyError: 


The VCF file needs to have columns #CHROM, POS … etc. The configuration script checks each VCF file for proper columns and will print a warning if any are missing or wrong. However, it does not halt execution and will include any malformed column VCF files and attempt to process them regardless. The main script may finish execution with the malformed VCF, but the output may be perturbed or useless.

  File "mucor/inputs.py", in parse_VarScan
    position = int(row[fieldId['POS']])
KeyError: 'POS'


Users may not supply vcf files that have inconsistent ‘effect’ and/or ‘functional consequence’ for the same mutation. Presumably, if the mutations has the same location and reference and alternate alleles, the effect and functional consequence must also be the same. This will cause a problem when collapsing mutations in the variant_details output type(s).
