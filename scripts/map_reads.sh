#!/bin/bash
#SBATCH -J amplicon_seq #job name
#SBATCH --time=0-1:30:00 
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32g
#SBATCH --output=genome_index.%j.out
#SBATCH --error=genome_index.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jacob.dayton@tufts.edu

mkdir demult
mkdir demult_merge
mkdir mapped
mkdir vcf

module load anaconda/2021.05 #or module load anaconda/2021.11
source activate ampliconez

current_dir=$(pwd)
LIBRARY_NAME=$(basename "$current_dir" | tr '[:lower:]' '[:upper:]')
echo "LIBRARY_NAME=$LIBRARY_NAME"

source activate /cluster/tufts/bio/tools/conda_envs/trim_galore/
trim_galore --paired --length 36 --dont_gzip -o trim ./*.fastq.gz --no_report_file

source deactivate

module load cutadapt/3.7

cutadapt -e 0.1 --no-indels \
    -g file:barcodes_fwd.fasta \
    -G file:barcodes_rev.fasta \
    -o ./demult/{name1}-{name2}.1.fastq.gz -p ./demult/{name1}-{name2}.2.fastq.gz \
    ./trim/*_R1*.fq ./trim/*_R2*.fq

rm -f ./demult/*unknown*.gz

for R1_file in ./demult/*.1.fastq.gz; do
    
    # Generate the corresponding R2 file name
    R2_file="${R1_file/.1.fastq.gz/.2.fastq.gz}"

    # Extract the sample name from the R1 file name
    sample=$(basename "$R1_file" .1.fastq.gz)

    #Merge paired-end reads using NGmerge
    NGmerge -1 $R1_file -2 $R2_file -o ./demult_merge/${LIBRARY_NAME}_${sample}.merged.fastq -v -n 2 -d -m 50
    
    #Filter out sequences < 110 bases
    zcat ./demult_merge/${LIBRARY_NAME}_${sample}.merged.fastq | seqtk seq -L 110 > ./demult_merge/${LIBRARY_NAME}_${sample}.selected.fastq

    echo "Reads corresponding to $sample have been merged"
done

#make hisat2 genome index in ampliconez parent directory
hisat2-build /cluster/tufts/dopmanlab/Jacob/ampliconez/reference_ampliconez3.fasta genome_index

for R1_file in ./demult_merge/*.selected.fastq; do

    # Extract the sample name from the R1 file name
    sample=$(basename "$R1_file" .selected.fastq)
    echo $sample
    hisat2 -x /cluster/tufts/dopmanlab/Jacob/ampliconez/genome_index -U ${R1_file} \
    -k 1 \
    --score-min L,0,-0.6 \
    --trim5 25 --trim3 25 \
    --no-spliced-alignment \
    --rg-id $sample --rg SM:$sample  --rg LB:lib1 --rg PU:unit1 --rg PL:illumina | \
    samtools view -bS - | \
    samtools sort -o ./mapped/${sample}.sorted.bam
    samtools index ./mapped/${sample}.sorted.bam

    echo "Alignment completed for $sample"
done

#separately run haplotype caller for males, females, and unknown sexed individuals

#haplogen_males.sh
#haplogen_females.sh
#haplogen_unknown.sh