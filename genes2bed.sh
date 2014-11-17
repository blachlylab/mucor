reference=${1} && genes=${*:2}
for gene in $genes;
do
    fgrep '"'$gene'"' $reference | awk -F'\t' '$3=="gene" {print $1"\t"$4"\t"$5}' | sort |uniq
done > out.bed