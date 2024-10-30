
library(data.table)
library(Matrix)

#Analysis of RNAseq data
#Generate read/count matrix from StarSolo output
meta <- read.table("./data/brbseq_combined.csv", header=TRUE,
                   sep = ',')
meta <- meta[,2:ncol(meta)]

##########################
#reorder matrix counts for map_d3_v2
matrix_dir <- "./data/map_d3_v2_norrna/Solo.out/Gene/raw/"
f <- file(paste0(matrix_dir, "umiDedup-NoDedup.mtx"), "r")
counts_d3 <- as.data.frame(as.matrix(readMM(f)))
close(f)

feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,
                      stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE,
                         stringsAsFactors = FALSE, data.table = F)
colnames(counts_d3) <- barcode.names$V1
rownames(counts_d3) <- feature.names$V1
rownames(meta) <- meta$code

#verify that all row labels in meta are columns in data
all(colnames(counts_d3) == meta$barcode)
idx = match(meta$barcode, colnames(counts_d3))
counts_d3 = counts_d3[,idx]
all(colnames(counts_d3) == meta$barcode)
colnames(counts_d3) <- meta$code
all(colnames(counts_d3) == meta$code)

fwrite(counts_d3, file = "./data/umi.counts_d3_v2.txt", sep = "\t", quote = F, row.names = T,
       col.names = T)
       
#reorder matrix counts for map_d2_v2
matrix_dir <- "./data/map_d2_v2_norrna/Solo.out/Gene/raw/"

f <- file(paste0(matrix_dir, "umiDedup-NoDedup.mtx"), "r")
counts_d2 <- as.data.frame(as.matrix(readMM(f)))
close(f)

feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,
                      stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE,
                      stringsAsFactors = FALSE, data.table = F)
colnames(counts_d2) <- barcode.names$V1
rownames(counts_d2) <- feature.names$V1
rownames(meta) <- meta$code

#verify that all row labels in meta are columns in data
all(colnames(counts_d2) == meta$barcode)
idx = match(meta$barcode, colnames(counts_d2))
counts_d2 = counts_d2[,idx]
all(colnames(counts_d2) == meta$barcode)
colnames(counts_d2) <- meta$code
all(colnames(counts_d2) == meta$code)

fwrite(counts_d2, file = "./data/umi.counts_d2_v2.txt", sep = "\t", quote = F, row.names = T,
       col.names = T)

##########################
## Gene-level differential expression analysis using BRBseq
library(RColorBrewer)
library(pheatmap)
library(readr)
library(edgeR)
library(vsn) #
library(dplyr)
library(ggrepel) 
library(data.table)
library(limma)
library(circacompare)
library(tidyverse)
library(ggpubr)

###################################################
counts_d2 <- read.table("./data/umi.counts_d2_v2.txt",
                        header=TRUE, row.names = 1)
counts_d3 <- read.table("./data/umi.counts_d3_v2.txt",
                        header=TRUE, row.names = 1)

#verify that all dimensions match between both counts files
all(rownames(counts_d2) == rownames(counts_d3))
all(colnames(counts_d2) == colnames(counts_d3))

#combine technical replicates
counts = cbind(counts_d2 + counts_d3)

#verify that colsums are the same
all(colSums(counts_d2)+colSums(counts_d3) == colSums(counts))
rm(counts_d2)
rm(counts_d3)
#verify that all row labels in meta are columns in data
all(colnames(counts) == rownames(meta))

#set variables as factors
meta$treatment = as.factor(meta$treatment)
meta$bkgrnd = as.factor(meta$bkgrnd)
meta$rep = as.factor(meta$rep)
meta$genotype = as.factor(meta$genotype)
meta$application = as.factor(meta$application)

all(colnames(counts) == rownames(meta))

# Compute depth per sample
meta$depth <- with(meta, colSums(counts))

table(rowSums(counts) != 0)

ggplot(data = meta %>% filter(treatment == '12L:12D', bkgrnd != 'bz', genotype == 'wildtype'),
       aes(x = treatment, y = depth/1000000, color = bkgrnd)) +
  geom_boxplot() + geom_point(position=position_jitterdodge(dodge.width=0.8)) +
  scale_y_continuous(breaks = seq(0,30, by=4), name = "depth (million reads)") +
  theme_bw() + geom_hline(yintercept=0.2, linetype = 'dashed')


# Filter samples that did not aggregate at least 0.5M reads
outliers <- rownames(subset(meta, depth < 0.5E6))
counts_sub <- counts[,-which(colnames(counts) %in% outliers)]

meta_sub <- meta[-which(colnames(counts) %in% outliers), ]

all(colnames(counts_sub) == rownames(meta_sub))

#set time as a factor: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments
meta_sub$azt = as.factor(meta_sub$azt) 


#filter to retain only wildtype 12:12 samples
meta_sub = meta_sub %>% filter(treatment == '12L:12D', genotype == 'wildtype', application == 'none',
                               bkgrnd != 'bz') %>% droplevels()
idx_sub = match(rownames(meta_sub), colnames(counts_sub))
counts_sub = counts_sub[, idx_sub]


# Filter low expressed genes
table(rowSums(counts_sub) != 0)
counts_sub <- counts_sub[rowSums(edgeR::cpm(counts_sub) >= 0.5) >= ncol(counts_sub)*0.2,] # At least 0.25 count per million in 20% samples

table(rowSums((counts_sub))==0) #any genes that are not expressed?
nrow(counts_sub)
with(meta_sub, table(bkgrnd, azt))

with(meta_sub, table(depth >= 1E6, azt, bkgrnd))
with(meta_sub, table(depth >= 0.5E6, azt, bkgrnd))

meta_sub %>% group_by(bkgrnd, treatment) %>%
  summarize(reads.mu = mean(depth),
            reads.med = median(depth),
            reads.max = max(depth),
            reads.min = min(depth),
            reads.lb = t.test(depth)$conf.int[1],
            reads.ub = t.test(depth)$conf.int[2])

#filter to remove D replicates
meta_sub = meta_sub %>% filter(code != "B1B", code != "B8B",
                               code != "F2D", code != "F3D",
                               code != "F4D", code != "F5D",
                               code != "F6D", code != "F7D",
                               code != "F8D")

idx_sub = match(rownames(meta_sub), colnames(counts_sub))
counts_sub = counts_sub[, idx_sub]

all(colnames(counts_sub) == rownames(meta_sub))

meta_sub = droplevels(meta_sub)


meta_sub %>% group_by(bkgrnd, azt) %>%
  summarize(depth.mu = mean(depth))

##################################
#Voom transforcounts_sub#Voom transformation and calculation of variance weights
#1. Counts are transformed to log2 counts per million reads (CPM), where “per million reads” is defined based on the normalization factors we calculated earlier
#2. A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
#3. A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
#4. The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.

mm <- with(meta_sub, model.matrix(~0+bkgrnd*azt))
#Transform count data to log2-counts per million (logCPM)
y_voom <- voom(counts_sub,mm,plot=T)

colnames(mm) = make.names(colnames(mm))

y_voom$E[1:5,1:2] #numeric matrix of normalized expression values on log2 scale

summary(prcomp(t(y_voom$E), center = T, scale. = T))
data.pca <- data.frame(prcomp(t(y_voom$E), center = T, scale. = T)$x)
data.pca <- cbind(data.pca, meta_sub)

# Plot
pca = ggplot(data.pca, aes(x = PC1, y = PC2, color = bkgrnd)) +
  geom_point(size=3, alpha = 0.75) +
  geom_text_repel(aes(label=azt), color='black', size = 5/14*8) +
  theme_bw()  +
  theme(text = element_text(size = 14, family = "sans", color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'top') +
  scale_x_continuous(breaks = seq(-30, 30, by=30), name = 'PC1: 37.6% variance',
                     expand=c(0,0), limits=c(-35,35)) +
  scale_y_continuous(breaks = seq(-30, 30, by=30), name = 'PC2: 34.3% variance',
                     expand=c(0,0), limits=c(-35,35))+
  scale_colour_manual(values = c('b_nd' = "#00BFC4", 'fh_long' = "#F8766D"),
                      labels = c('fh_long' = "Univoltine",
                                 'b_nd' = "Bivoltine"))

pca
ggsave("./plots/pca_expression_12L12D.png", plot = pca, width = 6, height = 4, dpi = 800, bg = 'white')

####################################################################################
#identify rhythmically expressed genes, genes that are differentially rhythmic, and genes that are differentially expressed
#using limorhyde

library('data.table')
library('foreach')
library('ggplot2')
library('limma')
library('limorhyde')
library('qs')

#specify zeitgeber per and q-value cutoffs for rhythmic and differentially rhythmic genes
per = 24
qvalRhyCutoff = 0.10
qvalDrCutoff = 0.1

#Convert data into a format that limorhyde base functions prefer :)

#Overwrite y to contain only the limma-voom normalized gene counts
y = y_voom$E
class(y) #"matrix" "array"
metadata = meta_sub

#Change "azt" to "time" and "bkgrnd" to "cond"
metadata$azt = as.numeric(as.character(metadata$azt))
metadata$azt
colnames(metadata)[11] = 'cond'
colnames(metadata)[11]
colnames(metadata)[13] = 'time'
colnames(metadata)

#Convert metadate "data.frame" to a "data.table" "data.frame"
metadata = as.data.table(metadata)

#Calculate time_cos and time_sin for each timepoint and append to the metadata data.table
metadata = cbind(metadata, limorhyde(metadata$time, 'time_'))

#Identify rhythmic genes: Calculate q-value of rhythmicity for each gene using that gene's p-value for each condition and adjusting for multiple testing
rhyLimma = foreach(condNow = unique(metadata$cond), .combine = rbind) %do% {
  design = model.matrix(~ time_cos + time_sin, data = metadata[cond == condNow])
  fit = lmFit(y[, metadata$cond == condNow], design)
  fit = eBayes(fit, trend = TRUE)
  rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
  setnames(rhyNow, 'rn', 'gene_id')
  rhyNow[, cond := condNow]}

rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
setorderv(rhyLimmaSummary, 'adj.P.Val')

nrow(rhyLimmaSummary) #q.rhythmic <0.15, Singer & Hughey 2019
mean(rhyLimmaSummary$adj.P.Val <= qvalRhyCutoff) # ~13.6% of genes were rhythmic
mean(rhyLimmaSummary$adj.P.Val <= 0.10) # ~13.6% of genes were rhythmic
mean(rhyLimmaSummary$adj.P.Val <= 0.05) # ~13.6% of genes were rhythmic

candidate_genes = read.csv("gene_ids.csv")
nrow(candidate_genes)

rhyLimmaSummary$rank = rownames(rhyLimmaSummary)

rhyLimmaSummary %>% filter(gene_id %in% candidate_genes$gene_id) %>%
  left_join(candidate_genes, by="gene_id")

rhyLimmaSummary_cg = rhyLimmaSummary %>% filter(gene_id %in% candidate_genes$gene_id) %>%
  left_join(candidate_genes, by="gene_id")

rhyLimmaSummary_cg$isRhy = rhyLimmaSummary_cg$adj.P.Val <= qvalRhyCutoff

rhyLimmaSummary_cg

######
#Identify differentially rhythmic genes, based on interactions between cond*time
#We pass all genes to limma (Empirical Bayes is best with many genes),
#but keep results only for rhythmic genes and adjust for multiple testing
design = model.matrix(~ cond * (time_cos + time_sin), data = metadata)

fit = lmFit(y, design)
fit = eBayes(fit, trend = TRUE)
drLimma = data.table(topTable(fit, coef = 5:6, number = Inf), keep.rownames = TRUE)
setnames(drLimma, 'rn', 'gene_id')

#subset only rhythmic genes (prior to check)
drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= qvalRhyCutoff]$gene_id]
drLimma[, adj.P.Val := round(p.adjust(P.Value, method = 'BH'),3)]
setorderv(drLimma, 'adj.P.Val')

drLimma$isDR <- (drLimma$adj.P.Val <= qvalDrCutoff)
table(drLimma$adj.P.Val <= qvalDrCutoff)

drLimma$rank = rownames(drLimma)

#q <= 0.1 is differentially rhythmic
drLimma_cg = drLimma %>% filter(gene_id %in% candidate_genes$gene_id) %>%
  left_join(candidate_genes, by="gene_id")


######
#Identify differentially expressed genes (based on interaction for cond with no interaction terms)
#Keep results only for non-differentially rhythmic genes

design = model.matrix(~ cond + time_cos + time_sin, data = metadata)

fit = lmFit(y, design)
fit = eBayes(fit, trend = TRUE)
deLimma = data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = TRUE)
setnames(deLimma, 'rn', 'gene_id')

nrow(deLimma)

deLimma = deLimma[!(gene_id %in% drLimma[adj.P.Val <= qvalDrCutoff]$gene_id)]
nrow(deLimma)
deLimma[, adj.P.Val := round(p.adjust(P.Value, method = 'BH'),3)]
setorderv(deLimma, 'adj.P.Val')

deLimma$rank = rownames(deLimma)

#differentially expressed genes
deLimma %>% filter(gene_id %in% candidate_genes$gene_id) %>%
  left_join(candidate_genes, by="gene_id")
#DEG <= 0.05

DEGs = deLimma[which(deLimma$adj.P.Val<0.01 & (deLimma$logFC > 0.5 | deLimma$logFC < -0.5)),]
nrow(DEGs)
View(DEGs)

DEGs_cg = DEGs %>% filter(gene_id %in% candidate_genes$gene_id) %>%
  left_join(candidate_genes, by="gene_id")


deLimma_cg = deLimma %>% filter(gene_id %in% candidate_genes$gene_id) %>%
  left_join(candidate_genes, by="gene_id")

deLimma_cg$rank = as.numeric(deLimma_cg$rank)

####################################################################################
cg_plots = list() #store plot for each candidate gene
cg_anova = list() #dataframe with a col of candidate_gene ids and col of 
cg_circacompare = list()
cg_names = vector()
cg_emms = list()

rhyLimmaSummary_cg$isRhy = TRUE
library(emmeans)
for (i in 1:nrow(rhyLimmaSummary_cg)) {
  gene_id = rhyLimmaSummary_cg[with(rhyLimmaSummary_cg, isRhy==TRUE)]$gene_id[i]
  gene_id_tab = data.frame(count = y_voom$E[which(rownames(y_voom$E) == gene_id),],
                           sample = rownames(meta_sub),
                           bkgrnd = meta_sub$bkgrnd,
                           azt = meta_sub$azt,
                           rep = meta_sub$rep)
  
  gene_id_tab$azt = as.factor(gene_id_tab$azt)
  
  print(max(gene_id_tab$count))
  
  m1 = lm(count ~ bkgrnd*azt,
          data = gene_id_tab)
  
  cg_anova[[i]] = car::Anova(m1)
  cg_emms[[i]] = pairs(emmeans(m1, ~ bkgrnd | azt,type = "response"), adjust = 'BH')
  
  gene_id_tab$azt = as.numeric(as.character(gene_id_tab$azt))
  
  gene_id_tab2 = gene_id_tab
  colnames(gene_id_tab2) = c("measure", "sample", "group", "time", "rep")
  
  title = rhyLimmaSummary_cg$gene_name[rhyLimmaSummary_cg$gene_id == gene_id]
  result2 <- circacompare(x = gene_id_tab2, col_time = "time", col_group = "group", col_outcome = "measure", per=24,
                          alpha_threshold=0.99)
  result2$summary[,2] = round(result2$summary[,2], 3)
  result2$plot <- result2$plot + theme_bw() +
    theme(text = element_text(size = 14, family = "sans", color = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = seq(2,23, by=3), name = 'Zeitgeber Time (hours)',
                       expand=c(0,0), limits=c(NA,25),
                       labels = c("14","17","20","23",
                                  "2","5","8", "11")) +
    scale_colour_manual(values = c('b_nd' = "#00BFC4", 'fh_long' = "#F8766D"),
                        labels = c('fh_long' = "Univoltine",
                                   'b_nd' = "Bivoltine")) +
    annotate("rect", xmin=0, xmax=12, ymin=-3, ymax=9, alpha=0.2, fill="gray") +
    scale_y_continuous(seq(-2,8,by=2), name = expression("Expression"*"(log"[2]*"CPM)"),
                       expand=c(0,0), limits=c(-3,9)) +
    ggtitle(as.expression(bquote(italic(.(rhyLimmaSummary_cg$gene_name[rhyLimmaSummary_cg$gene_id == gene_id])))))
  
  result2$plot$layers[[3]]$aes_params$size = 2
  result2$plot$layers[[3]]$aes_params$alpha = 0.5
  cg_plots[[i]] = result2$plot
  cg_circacompare[[i]] = result2$summary
  cg_names[i] = rhyLimmaSummary_cg$gene_name[i]
  i = i+1
}

names(cg_plots) = unlist(cg_names)
names(cg_circacompare) = unlist(cg_names)
names(cg_anova) = unlist(cg_names)
names(cg_emms) = unlist(cg_names)


#clk
cg_plots$clk = cg_plots$clk + scale_y_continuous(seq(-3,3,by=3),
                                                     name = expression("Expression"*"(log"[2]*"CPM)"),
                                                     expand=c(0,0),
                                                     limits=c(-4.5,4.5)) +
  annotate("rect", xmin=0, xmax=12, ymin=-4.5, ymax=4.5, alpha=0.2, fill="gray")

cg_plots$clk$layers[[1]]$aes_params$linetype = 'longdash'


#bmal1  
cg_plots$bmal1 = cg_plots$bmal1 + scale_y_continuous(seq(0,6,by=3),
                                                     name = expression("Expression"*"(log"[2]*"CPM)"),
                                                     expand=c(0,0),
                                                     limits=c(-1.5,7.5)) +
  annotate("rect", xmin=0, xmax=12, ymin=-1.5, ymax=7.5, alpha=0.2, fill="gray")

cg_plots$bmal1$layers[[2]]$aes_params$linetype = 'longdash'


#cyc
cg_plots$cyc = cg_plots$cyc + scale_y_continuous(seq(0,6,by=3),
                                                     name = expression("Expression"*"(log"[2]*"CPM)"),
                                                     expand=c(0,0),
                                                     limits=c(-1.5,7.5)) +
  annotate("rect", xmin=0, xmax=12, ymin=-1.5, ymax=7.5, alpha=0.2, fill="gray")

cg_plots$cyc$layers[[1]]$aes_params$linetype = 'longdash'


#cry2
cg_plots$cry2 = cg_plots$cry2 + scale_y_continuous(seq(2,6,by=2),
                                                   name = expression("Expression"*"(log"[2]*"CPM)"),
                                                   expand=c(0,0),
                                                   limits=c(1,7)) +
  annotate("rect", xmin=0, xmax=12, ymin=1, ymax=7, alpha=0.2, fill="gray")

cg_plots$cry2 = cg_plots$cry2 +
  annotate("segment",
           x = cg_circacompare$cry2[11,2],
           xend = cg_circacompare$cry2[11,2],
           y = 1,
           yend = sum(cg_circacompare$cry2[3,2]+cg_circacompare$cry2[7,2]),
           linetype="dashed", color = "black", linewidth = 0.5) +
  annotate("segment",
           x = cg_circacompare$cry2[12,2],
           xend = cg_circacompare$cry2[12,2],
           y = 1,
           yend = sum(cg_circacompare$cry2[4,2]+cg_circacompare$cry2[8,2]),
           linetype="dashed", color = "black", linewidth = 0.5)


#per
cg_plots$per = cg_plots$per + scale_y_continuous(seq(2,6,by=2),
                                                       name = expression("Expression"*"(log"[2]*"CPM)"),
                                                       expand=c(0,0),
                                                       limits=c(1,7)) +
  annotate("rect", xmin=0, xmax=12, ymin=1, ymax=7, alpha=0.2, fill="gray")


cg_plots$per = cg_plots$per +
  annotate("segment",
           x = cg_circacompare$per[11,2],
           xend = cg_circacompare$per[11,2],
           y = 1,
           yend = sum(cg_circacompare$per[3,2]+cg_circacompare$per[7,2]),
           linetype="dashed", color = "black", linewidth = 0.5) +
  annotate("segment",
           x = cg_circacompare$per[12,2],
           xend = cg_circacompare$per[12,2],
           y = 1,
           yend = sum(cg_circacompare$per[4,2]+cg_circacompare$per[8,2]),
           linetype="dashed", color = "black", linewidth = 0.5)

#cwo
cg_plots$cwo$layers[[5]] <- NULL #delete rectangle
cg_plots$cwo = cg_plots$cwo + scale_y_continuous(seq(-5,5,by=5),
                                                 name = expression("Expression"*"(log"[2]*"CPM)"),
                                                 expand=c(0,0),
                                                 limits=c(-7.5,7.5)) +
  annotate("rect", xmin=0, xmax=12, ymin=-7.5, ymax=7.5, alpha=0.2, fill="gray")

#make line for fh_long dashed because this group is not rhythmic
cg_plots$cwo$layers[[2]]$aes_params$linetype = 'longdash'

#tim
cg_plots$tim = cg_plots$tim + scale_y_continuous(seq(4,8,by=2),
                                                           name = expression("Expression"*"(log"[2]*"CPM)"),
                                                           expand=c(0,0),
                                                           limits=c(3,9)) +
  annotate("rect", xmin=0, xmax=12, ymin=3, ymax=9, alpha=0.2, fill="gray")

#make line for fh_long dashed because this group is not rhythmic
cg_plots$tim$layers[[2]]$aes_params$linetype = 'longdash'

#pdp1-e
cg_plots$`pdp1-e` = cg_plots$`pdp1-e`  + scale_y_continuous(seq(4,8,by=2),
                                                            name = expression("Expression"*"(log"[2]*"CPM)"),
                                                            expand=c(0,0),
                                                            limits=c(3,9)) +
  annotate("rect", xmin=0, xmax=12, ymin=3, ymax=9, alpha=0.2, fill="gray")

#make line for fh_long dashed because this group is not rhythmic
cg_plots$`pdp1-e`$layers[[2]]$aes_params$linetype = 'longdash'


#vri
cg_plots$vri = cg_plots$vri  + scale_y_continuous(seq(2,6,by=2),
                                                        name = expression("Expression"*"(log"[2]*"CPM)"),
                                                        expand=c(0,0),
                                                        limits=c(1,7)) +
  annotate("rect", xmin=0, xmax=12, ymin=1, ymax=7, alpha=0.2, fill="gray")

cg_plots$vri = cg_plots$vri +
  annotate("segment",
           x = cg_circacompare$vri[11,2],
           xend = cg_circacompare$vri[11,2],
           y = 1,
           yend = sum(cg_circacompare$vri[3,2]+cg_circacompare$vri[7,2]),
           linetype="dashed", color = "black", linewidth = 0.5) +
  annotate("segment",
           x = cg_circacompare$vri[12,2],
           xend = cg_circacompare$vri[12,2],
           y = 1,
           yend = sum(cg_circacompare$vri[4,2]+cg_circacompare$vri[8,2]),
           linetype="dashed", color = "black", linewidth = 0.5) 



cg_plots$vri = cg_plots$vri + annotate("text", x = 6, y = 2, label = expression(Phi~"= -1.21 hours"), size = 3.5)
cg_plots$per = cg_plots$per + annotate("text", x = 6, y = 6.5, label = expression(Phi~"= -1.04 hours"), size = 3.5)
cg_plots$cry2 = cg_plots$cry2 + annotate("text", x = 6, y = 6.5, label = expression(Phi~"= -1.38 hours"), size = 3.5)



final_plot = ggarrange(cg_plots$cry2 + rremove("xlab") + rremove("ylab"),
                       cg_plots$per + rremove("xlab") + rremove("ylab"),
                       cg_plots$vri + rremove("xlab") + rremove("ylab"),
                       cg_plots$tim + rremove("xlab") + rremove("ylab"),
                       cg_plots$`pdp1-e` + rremove("xlab") + rremove("ylab"),
                       cg_plots$cwo + rremove("xlab") + rremove("ylab"),
                       ncol=3, nrow=2,
                       common.legend = TRUE, align = 'h')

library(grid)

final_plot = annotate_figure(final_plot, left = textGrob(expression("Expression"*"(log"[2]*"CPM)"), rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                             bottom = textGrob("Zeitgeber Time (hours)", gp = gpar(cex = 1.3)))

final_plot
#ggsave("clk_targets.png", plot = final_plot, width = 10, height = 6, dpi = 800, bg = 'white')


all_ccgs = ggarrange(cg_plots$clk + rremove("xlab") + rremove("ylab"),
                     cg_plots$bmal1 + rremove("xlab") + rremove("ylab"),
                     cg_plots$cyc + rremove("xlab") + rremove("ylab"),
                     cg_plots$cry2 + rremove("xlab") + rremove("ylab"),
                     cg_plots$per + rremove("xlab") + rremove("ylab"),
                     cg_plots$vri + rremove("xlab") + rremove("ylab"),
                     cg_plots$tim + rremove("xlab") + rremove("ylab"),
                     cg_plots$`pdp1-e` + rremove("xlab") + rremove("ylab"),
                     cg_plots$cwo + rremove("xlab") + rremove("ylab"),
                     ncol=3, nrow=3,
                     common.legend = TRUE, align = 'h')

all_ccgs = annotate_figure(all_ccgs, left = textGrob(expression("Expression"*"(log"[2]*"CPM)"), rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                           bottom = textGrob("Zeitgeber Time (hours)", gp = gpar(cex = 1.3)))

ggsave("./plots/figure2.png", plot = all_ccgs, width = 7.5, height = 6, dpi = 800, bg = 'white')

rhyLimmaSummary_cg$isRhy = rhyLimmaSummary_cg$adj.P.Val <= qvalRhyCutoff
rhyLimmaSummary_cg

cg_circacompare$cry2
cg_circacompare$per
cg_circacompare$vri
cg_circacompare$tim

sink(file = "cg_output_circa_anova_etc_STAR_v2.txt")
print("Summary of Rhythmic Genes")
rhyLimmaSummary_cg
print("")
print("Summary of Differentially Rhythmic Genes")
drLimma_cg
print("")
print("Summary of DEGs")
deLimma_cg
print("Top 15 DEGs and bottom 15 DEGs")
head(DEGs, 15)
tail(DEGs, 15)
print("Circacompare results")
cg_circacompare

print('results from read counts')

meta_sub %>% group_by(bkgrnd, treatment) %>%
  summarize(reads.mu = mean(depth),
            reads.med = median(depth),
            reads.max = max(depth),
            reads.min = min(depth),
            reads.lb = t.test(depth)$conf.int[1],
            reads.ub = t.test(depth)$conf.int[2])

sink()