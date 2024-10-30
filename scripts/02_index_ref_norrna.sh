#!/bin/bash
#SBATCH -J genome_index #job name
#SBATCH --time=0-16:00:00 
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=96g
#SBATCH --output=genome_index.%j.out
#SBATCH --error=genome_index.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jacob.dayton@tufts.edu

module load star/2.7.11b
STAR --runMode genomeGenerate --genomeDir /cluster/tufts/dopmanlab/Jacob/reference/no_rrna_star/ \
--genomeFastaFiles /cluster/tufts/dopmanlab/Jacob/reference/no_rrna_star/GCF_963855985.1_ilOstNubi1.1_genomic.fa --sjdbGTFfile \
/cluster/tufts/dopmanlab/Jacob/reference/no_rrna_star/GCF_genomic_no_rRNA.gtf --runThreadN 12 --genomeSAindexNbases 13
