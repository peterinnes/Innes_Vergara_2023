assemblies=(C24.phase0.fa cannatonic_GCA_001865755.1_ASM186575v1_genomic.fna GCF_900626175.2_cs10_genomic.fna C_sativa_isolate_JL_genome_only_chrs.fasta GCA_000230575.5_ASM23057v5_genomic.fna GCA_003417725.2_ASM341772v2_genomic.fna GCA_003660325.2_Oct15_3.7Mb_N50_Jamaican_Lion_Assembly_genomic.fna steepAsm-Quiver.fasta_unwrapped.fa USO31.phase0.fa GCA_000830395.1_hl_cordifolius_version_1.0.fasta_genomic.fna GCA_000831365.1_hl_lupulus_version_1.0.fasta_genomic.fna)
assembly_names=(C24 Cannatonic cs10 isolate_JL Purple_Kush Finola Jamaican_Lion PBBK USO31 Hl_cordifolius Hl_lupulus)

queries=(exons_LOC115713215_CsPT1 exons_LOC115713185_CsPT4 exons_LOC115713205_CsPT7)

rundir=~/CsPTx_genomics/data/getfasta_results/
outdir=~/CsPTx_genomics/data/blast_results/exon_blasts/

#IFS=$'\n' read -d '\n' -r -a array < assembly_paths.txt
#IFS=$'\n' read -d '\n' -r -a array2 < assembly_paths.txt

cd $rundir
for i in "${!assemblies[@]}"; do
    echo ${assemblies[i]} 

    if [[ -e ~/CsPTx_genomics/data/assemblies/${assemblies[i]}.nhr ]]; then
        echo "blast index exists"
    else
        echo "index does not exist... indexing now..."
        makeblastdb -in ~/CsPTx_genomics/data/assemblies/${assemblies[i]} -dbtype nucl
    fi
    
    for j in "${!queries[@]}"; do
        echo "running blastn:" ${queries[j]}
        blastn -query ${queries[j]}.cs10gff.fasta -db ~/CsPTx_genomics/data/assemblies/${assemblies[i]} -outfmt 6 | sort -k 9 -n > ${outdir}${assembly_names[i]}_${queries[j]}_2022-05-26.blast
    done
done
