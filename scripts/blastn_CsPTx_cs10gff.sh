ref=~/CsPTx_genomics/data/assemblies/GCF_900626175.2_cs10_cds_from_genomic.fna
queries=(CsPT1 CsPT4 CsPT7)

rundir=~/CsPTx_genomics/
outdir=~/CsPTx_genomics/data/blast_results/

#IFS=$'\n' read -d '\n' -r -a array < assembly_paths.txt
#IFS=$'\n' read -d '\n' -r -a array2 < assembly_paths.txt

cd $rundir

if [[ -e $ref.nhr ]]; then
    echo "blast index exists"
else
    echo "index does not exist... indexing now..."
    makeblastdb -in $ref -dbtype nucl
fi
    
for j in "${!queries[@]}"; do
    echo "running blastn:" ${queries[j]}
    blastn -query ../${queries[j]}.fasta -db $ref -outfmt 7 | sort -k 7 -n > ${outdir}cs10gff_${queries[j]}_2022-3-27.blast
done
