#!/bin/sh 
echo downloading data files
wget -i files_to_get.txt
echo decompressing annotation
gunzip gencode.v19.annotation.gtf.gz
echo finished
echo reformatting and indexing 00-common_all.vcf.gz
zcat 00-common_all.vcf.gz | bgzip > 00-common_all.vcf.bgz && rm 00-common_all.vcf.gz
mv 00-common_all.vcf.bgz 00-common_all.vcf.gz
tabix 00-common_all.vcf.gz
echo finished
echo setting up symbolic links for example runs
for i in 00-common_all.vcf.gz 00-common_all.vcf.gz.tbi gencode.v19.annotation.gtf;
do
    for j in example_workflow_1/projA example_workflow_1/projB example_workflow_2 example_workflow_3;
    do
	ln -s $(readlink -f ${i}) ${j}/
    done
done
