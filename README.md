# Period manuscript

**Description:** Scripts to reproduce the analyses reported in Dayton & Dopman ([2024](https://www.biorxiv.org/content/10.1101/2024.11.02.621642v1)), The *period* gene regulates daily and seasonal timing in *Ostrinia nubilalis*.

**Purpose:** Prior research in the European corn borer (*Ostrinia nubilalis*) associated ecotype differences in seasonal post-diapause development (PDD) time, a key factor determining voltinism (number of generations/growing season), with genetic variation in the circadian clock gene *period*. Here, I integrated behavioral assays, transcriptome profiling, and targeted genetic/pharmacological manipulations to investigate *per*’s contribute to both daily and seasonal biological rhythms.

**Graphic Summary of Key Results:**
<p align="center">
  <img src="https://github.com/user-attachments/assets/386cc7fb-a045-4784-8357-d2b6120d4e2a" alt="graphic"/>
</p>

***per* transcript abundance alters daily and seasonal rhythms in O. nubilalis.** (A) Variation in the phase (Φ) or amplitude of *per* mRNA abundance, and predicted PER protein accumulation, across a 24-hr day is associated with differences in daily behavior and seasonal development of univoltine (red) and bivoltine (blue) individuals. (B) In a simplified model of the Lepidopteran circadian clock network, we hypothesize that earlier vs. later repression of CLK:BMAL1 in the evening (scotophase) pleiotropically alters both daily and seasonal rhythms of wild-type individuals. (C) Experimental reductions in *per* gene dosage (blue) delay the phase of daily behavior, decrease photoperiodic diapause incidence, and shorten post-diapause development (PDD).

<ins>Workflow:</ins>
  1) Processing Amplicon-seq data to genotype mutant vs. wildtype individuals (GATK pipeline; to account for different zygosity of Z chr, call variants separately m v. f)
     - map_reads.sh (input: .fastq, barcodes_fwd.fasta, barcodes_rev.fasta; hisat2)
     - haplogen_females.sh & haplogen_males.sh & haplogen_unknown.sh (input: male_seqs.txt, female_seqs.txt, unknown_seqs.txt; GATK)
     - **Output:** VCFs
  2) Processing RNA-seq data (multiplexed 3'-end RNA-seq with inline barcode) to identify cycling genes & DEG circadian clock genes
     - 01_remove_rrna_from_gff.sh (grep)
     - 02_index_ref_norrna.sh (STAR)
     - 03_map_reads_norrna.sh (STAR/STARsolo) 
     - analysis_brbseq.R (limma-voom & edgeR, limorhyde)
     - **Output:** Gene/Count Matrix & DEG analysis
  3) Phenotypic differences among wild-type individuals
     - analysis_ecotypes.Rmd
     - analysis_interruptions.Rmd
     - analysis_diapauseXphoto.Rmd
     - analysis_licl.Rmd
     - **Outputs:** GLMs & figs
  5) Phenotype differences among mutant vs. wild-type individuals
     - analysis_interruptions.Rmd
     - analysis_per.Rmd
     - **Outputs:** GLMs & figs
     
