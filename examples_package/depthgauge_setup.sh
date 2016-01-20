## optionally fetch 100+ GB of BAM files for depth gauge example
cd ./data/
echo "Fetching data ... "
wget -i ../bams_to_get.txt # this may need to be customized depending on network configuration 
echo "Fetching Complete!"
echo "Setting up example_workflow_5 folder"
[ ! -d ../example_workflow_5 ] && mkdir ../example_workflow_5
cd ../example_workflow_5
for fn in $(awk -F'/' '{print $10}' ../bams_to_get.txt);
do 
    [ ! -d ${fn%%.*} ] && mkdir ${fn%%.*}
    ln -s ../../data/$fn ./${fn%%.*}/${fn%%.*}${fn##*[0-9]}
    [ ! -e ./${fn%%.*}/${fn%%.*}.VarScan.snps.vcf ] && ln -s ../../data/${fn%%.*}.VarScan.snps.vcf ./${fn%%.*}/
done
for fn in example5.json hotspots.bed samples.txt; 
do 
    [ ! -e ./$fn ] && ln -s ../data/${fn} . 
done
echo "Depth Gauge example setup complete!"
