assembly_paths=(/data1/superstar/assembly_C24/Carmagnola_Phase0/C24.phase0.fa /data1/superstar/cannatonic_assembly/cannatonic_GCA_001865755.1_ASM186575v1_genomic.fna /data0/CannabisGenomics2021/assemblies/cs10/GCF_900626175.2_cs10_genomic.fna /home/superstar/Cannabis_synteny/C_sativa_isolate_JL_genome_only_chrs.fasta /data1/superstar/PK_assembly/GCA_000230575.5_ASM23057v5_genomic.fna /data1/superstar/finola_laverty_assembly/GCA_003417725.2_ASM341772v2_genomic.fna /data1/superstar/jl_assembly/GCA_003660325.2_Oct15_3.7Mb_N50_Jamaican_Lion_Assembly_genomic.fna /data1/superstar/assembly_PBBK/steepAsm-Quiver.fasta_unwrapped.fa /data1/superstar/assembly_USO31/USO31_Phase0/USO31.phase0.fa /data1/superstar/humulus_assemblies/GCA_000830395.1_hl_cordifolius_version_1.0.fasta_genomic.fna /data1/superstar/humulus_assemblies/GCA_000831365.1_hl_lupulus_version_1.0.fasta_genomic.fna)
 
assembly_names=(C24 Cannatonic cs10 isolate_JL Purple_Kush Finola Jamaican_Lion PBBK USO31 Hl_cordifolius Hl_lupulus)

queries=(CsPT1 CsPT2 CsPT3 CsPT4 CsPT5 CsPT6 CsPT7)

rundir=/home/hemp2021/student_folders/Everybody/petah_pahka/
outdir=/home/hemp2021/student_folders/Everybody/petah_pahka/blast_results/

#IFS=$'\n' read -d '\n' -r -a array < assembly_paths.txt
#IFS=$'\n' read -d '\n' -r -a array2 < assembly_paths.txt

cd $rundir
for i in "${!assembly_paths[@]}"; do
    echo ${assembly_paths[i]} 

    if [[ -e ${assembly_paths[i]}.nhr ]]; then
        echo "blast index exists"
    else
        echo "index does not exist... indexing now..."
        makeblastdb ${assembly_paths[i]}
    fi
    
    for j in "${!queries[@]}"; do
        echo "running blastn:" ${queries[j]}
        blastn -query ${queries[j]}.fasta -db ${assembly_paths[i]} -outfmt 6 | sort -k 7 -n > ${outdir}${assembly_names[i]}_${queries[j]}_2021-11-22.blast
    done
done
