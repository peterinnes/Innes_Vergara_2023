assembly_paths=(~/student_folders/assemblies/C24/C24_phase1_assembly.fasta ~/student_folders/assemblies/cannatonicGCA/cannatonic_GCA_001865755.1_ASM186575v1_genomic.fna ~/student_folders/assemblies/cs10/GCF_900626175.2_cs10_genomic.fna ~/student_folders/assemblies/CSativaIsolateJLGenome/C_sativa_isolate_JL_genome.fasta ~/student_folders/assemblies/GC_000230575/GCA_000230575.5_ASM23057v5_genomic.fna ~/student_folders/assemblies/GCA_000830395/GCA_000830395.1_hl_cordifolius_version_1.0.fasta_genomic.fna ~/student_folders/assemblies/GCA_000831365/GCA_000831365.1_hl_lupulus_version_1.0.fasta_genomic.fna ~/student_folders/assemblies/GCA_003417725/GCA_003417725.2_ASM341772v2_genomic.fna ~/student_folders/assemblies/JamaicanLion/GCA_003660325.2_Oct15_3.7Mb_N50_Jamaican_Lion_Assembly_genomic.fna ~/student_folders/assemblies/steepAsm/steepAsm-Quiver.fasta_unwrapped.fa ~/student_folders/assemblies/USO31_phase1_assembly/USO31_phase1_assembly.fasta)
 
assembly_names=(C24 Cannatonic cs10 C_sativa_isolate_JL Purple_Kush hl_cordifolius hl_lupulus Finola Jamaican_Lion PBBK USO31)

query=/data0/CannabisGenomics2021/Everybody/petah_pahka/geranyltransferase.fasta
outdir=/home/hemp2021/student_folders/Everybody/petah_pahka/blast_results/

#IFS=$'\n' read -d '\n' -r -a array < assembly_paths.txt
#IFS=$'\n' read -d '\n' -r -a array2 < assembly_paths.txt

for i in "${!assembly_paths[@]}"; do
    echo ${assembly_paths[i]} 

if [[ -e ${assembly_paths[i]}.nhr ]]; then
        echo "index exists"
    else
        echo "index does not exist... indexing now..."
        makeblastdb ${assembly_paths[i]}
    fi

    blastn -query $query -db ${assembly_paths[i]} -outfmt 6 | sort -k 7 -n > ${outdir}${assembly_names[i]}_CBGAs.blast

done
