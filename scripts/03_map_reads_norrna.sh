#!/bin/bash
#SBATCH -J map_d2_v2 #job name
#SBATCH --time=0-12:00:00 
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=128g
#SBATCH --output=map_d2_v2.%j.out
#SBATCH --error=map_d2_v2.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jacob.dayton@tufts.edu

module load star/2.7.11b

mappedDIR="map_d2_v2_norrna/"

mkdir $mappedDIR


STAR --runMode alignReads --outSAMmapqUnique 60 --runThreadN 16 \
--outSAMunmapped Within --soloStrand Forward --quantMode GeneCounts \
--outBAMsortingThreadN 8 --genomeDir /cluster/tufts/dopmanlab/Jacob/reference/no_rrna_star/ \
--soloType CB_UMI_Simple --soloCBstart 1 \
--outFilterScoreMinOverLread 0.25 \
--outFilterMatchNminOverLread 0.25 \
--limitBAMsortRAM 6303886190 \
--soloCBlen 14 --soloUMIstart 15 --soloUMIlen 14 --soloUMIdedup NoDedup 1MM_All \
--soloCellFilter None --soloCBwhitelist barcodes.txt --soloBarcodeReadLength 0 \
--soloFeatures Gene --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
--outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $mappedDIR --readFilesIn library2_dayton_2_S1_R2_001.fastq.gz library2_dayton_2_S1_R1_001.fastq.gz
