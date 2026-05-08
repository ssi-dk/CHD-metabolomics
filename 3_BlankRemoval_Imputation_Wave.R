# BlankRemoval_Imputation_Wave.R
# Author: Filip Ottosson
#
# This script performed filtering, imputation and normalization for the metabolomic characterization of 
# congenital heart disorders. The input file contains raw intensity values for all metabolites detected 
# in at least one of the samples in the study. This also includes metabolites that miss MS2 spectra. 
# In the only prior filtering step, redundant adduct ions have been excluded. 
# See the python script "2_filter_IIMN_adducts" for more details.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(ggpubr)
  library(cowplot)
  library(vegan)
  library(WaveICA2.0)
})

# ---------------------------------------------------------------------------
# Paths — edit these before running
# ---------------------------------------------------------------------------

feature_table_file <- "path/to/feature_table_IIMN_adduct_filtered.csv"
metadata_file      <- "path/to/metadata.csv"
output_dir         <- "path/to/output_directory"

# ---------------------------------------------------------------------------
# Import Data
# ---------------------------------------------------------------------------

ft        <- read.csv(feature_table_file)
met_names <- ft$row.ID
ft        <- ft[, 14:(ncol(ft) - 1)]   # select columns with metabolite data
dim(ft)

md <- read.csv(metadata_file)

# ---------------------------------------------------------------------------
# Calculate Sparsity
# ---------------------------------------------------------------------------

sum(ft == 0) / (ncol(ft) * nrow(ft))

sparsity <- cbind(colSums(ft == 0), colSums(ft != 0), colSums(ft != 0) / nrow(ft))
colnames(sparsity) <- c("N missing", "N not missing", "% not missing")
sparsity

sparsity[grep("PO", rownames(sparsity)), ]

# ---------------------------------------------------------------------------
# Blank removal
# ---------------------------------------------------------------------------
# In this step metabolite measurements that are not sufficiently more intense than what is measured 
# in the blank samples are removed. First, the average blank intensities are calculated for each metabolite. 
# For each metabolite, all measurments that fall below five times the respective blank intensity is transformed to zero.

ft <- t(ft)
rownames(ft) <- gsub('.Peak.area', '', rownames(ft))
colnames(ft) <- met_names
class(ft) <- "numeric"

if (length(which(is.na(rowSums(ft) == T))) != 0) {
  ft <- ft[-which(is.na(rowSums(ft) == T)), ]
}

blankmeans <- colMeans(ft[grep('PB', rownames(ft)), ])  # mean paper blank intensity per metabolite

blankids <- grep('PB', rownames(ft))  # position of paper blanks
blankids

for (i in 1:ncol(ft)) {
  ft[-blankids, i][which(ft[-blankids, i] < 5 * blankmeans[i])] <- 0  # replace low intensity measurements with zero
}

if (length(which(colSums(ft) == 0)) != 0) {
  ft <- ft[, -which(colSums(ft[-blankids, ]) == 0)]  # remove empty columns
}

dim(ft)

sum(blankmeans != 0) / sum(blankmeans == 0)

hist(colSums(ft == 0))

# ---------------------------------------------------------------------------
# Filter missing features
# ---------------------------------------------------------------------------

# After blank filtering, metabolites will be divided into two datasets. ft_75 contains the metabolites with 
# at least 75 % presence (max 25 % sparsity). ft_25 contains metabolites with a sparsity between 25-75 %. 
# We can consider ft_75 to containt metabolites with a continuous distribution, while ft_25 can be treated as binary variables. 

ft_75 <- ft[, colSums(ft[grep("_S", rownames(ft)), ] == 0) < (nrow(ft[grep("_S", rownames(ft)), ]) * 0.25)]
ft_25 <- ft[, colSums(ft[grep("_S", rownames(ft)), ] == 0) < (nrow(ft[grep("_S", rownames(ft)), ]) * 0.75)]
ft_25 <- ft_25[, !colnames(ft_25) %in% colnames(ft_75)]

dim(ft)
dim(ft_75)
dim(ft_25)

ft_75 <- cbind.data.frame(gsub('[0-9]+', '', sapply(1:nrow(ft_75), function(n) { strsplit(rownames(ft_75), "_")[[n]][3] })),       
                          as.numeric(as.character(gsub('[^0-9.-]', '', sapply(1:nrow(ft_75), function(n) { strsplit(rownames(ft_75), "_")[[n]][3] })))),
                          sapply(1:nrow(ft_75), function(n) { strsplit(rownames(ft_75), "_")[[n]][2] }),
                          ft_75)
colnames(ft_75)[1:3] <- c("SampleType", "RunOrder", "Batch")

head(ft_75)

ft_25 <- cbind.data.frame(gsub('[0-9]+', '', sapply(1:nrow(ft_25), function(n) { strsplit(rownames(ft_25), "_")[[n]][3] })),
                          as.numeric(as.character(gsub('[^0-9.-]', '', sapply(1:nrow(ft_25), function(n) { strsplit(rownames(ft_75), "_")[[n]][3] })))),
                          sapply(1:nrow(ft_25), function(n) { strsplit(rownames(ft_25), "_")[[n]][2] }),
                          ft_25)
colnames(ft_25)[1:3] <- c("SampleType", "RunOrder", "Batch")

head(ft_25)

# ---------------------------------------------------------------------------
# PCA Overview
# ---------------------------------------------------------------------------

set.seed(1234)
Data  <- ft_75
start <- 4
pca_mod <- prcomp(Data[, start:ncol(Data)], center = TRUE, scale = TRUE)

PC1 <- pca_mod$x[, 1]
PC2 <- pca_mod$x[, 2]
PC3 <- pca_mod$x[, 3]
PC4 <- pca_mod$x[, 4]
PC5 <- pca_mod$x[, 5]
plot(pca_mod)

pcplotdata <- cbind.data.frame(Data$SampleType, Data$RunOrder, Data$Batch, PC1, PC2, PC3, PC4, PC5)
colnames(pcplotdata) <- c("SampleType", "RunOrder", "Batch", "PC1", "PC2", "PC3", "PC4", "PC5")

p_PC1_2 <- ggplot(pcplotdata, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point() + theme_classic() + scale_color_jama()

p_PC3_4 <- ggplot(pcplotdata, aes(x = PC3, y = PC4, color = SampleType)) +
  geom_point() + theme_classic() + scale_color_jama()

p_PC1_2_batch_samp <- ggplot(subset(pcplotdata, SampleType == "S" | SampleType == "PB"),
                              aes(x = PC1, y = PC2, color = Batch, shape = SampleType)) +
  geom_point() + theme_classic() + scale_color_jama()

p_PC1_2_batch <- ggplot(subset(pcplotdata, !SampleType == "S"),
                         aes(x = PC1, y = PC2, color = Batch, shape = SampleType)) +
  geom_point() + theme_classic() + scale_color_jama()

p_PC3_4_batch <- ggplot(subset(pcplotdata, !SampleType == "S"),
                         aes(x = PC3, y = PC4, color = Batch, shape = SampleType)) +
  geom_point() + theme_classic() + scale_color_jama()

ggsave(file.path(output_dir, "QC/PCA_SampleBatch.pdf"), p_PC1_2_batch_samp)
ggsave(file.path(output_dir, "QC/PCA_SampleType.pdf"), plot_grid(p_PC1_2, p_PC3_4, nrow = 2))
ggsave(file.path(output_dir, "QC/PCA_Batch.pdf"), plot_grid(p_PC1_2_batch, p_PC3_4_batch, nrow = 2))

p_PC1_2_batch_samp
plot_grid(p_PC1_2, p_PC3_4, nrow = 2)
plot_grid(p_PC1_2_batch, p_PC3_4_batch, nrow = 2)

# ---------------------------------------------------------------------------
# RSD in ECs and Pools
# ---------------------------------------------------------------------------
# Some metabolites will likely be quite noisy, i.e. the technical variation is intolerable. 
# We can filter these out by estimating the technical variation by looking at the stability 
# of the intensity measurements across repeated injections of external controls (EC) and plate pools (PO). 
# If the relative standard deviation is exceeding 30 % in both the ECs and in the POs, the metabolite 
# is filtered out. The Pool RSD is calculated per plate, since the pools on different plates are different 
# (pool of all samples in that plate). 

EC     <- subset(Data, SampleType == "EC")[, start:ncol(Data)]
RSD_EC <- sapply(1:ncol(EC), function(n) { sd(EC[, n]) / mean(EC[, n]) })
hist(RSD_EC)
sum(RSD_EC < 0.3, na.rm = TRUE)

PO   <- subset(Data, SampleType == "PO")[, start:ncol(Data)]
PO_1 <- PO[grep("CHD_1", rownames(PO)), ]

RSD_pool <- sapply(1:ncol(PO_1), function(n) { sd(PO_1[, n]) / mean(PO_1[, n]) })
hist(RSD_pool)
sum(RSD_pool < 0.3, na.rm = TRUE)

# ---------------------------------------------------------------------------
# Filter by RSD
# ---------------------------------------------------------------------------

RSD_pool[is.na(RSD_pool)] <- 1
keep_po <- names(PO[, RSD_pool < 0.3])

RSD_EC[is.na(RSD_EC)] <- 1
keep_ec <- names(EC[, RSD_EC < 0.3])

keep <- c(keep_po, keep_ec)
keep <- keep[!duplicated(keep)]
length(keep)

# ---------------------------------------------------------------------------
# Impute missing values
# ---------------------------------------------------------------------------
# Impute missing values by replacing with a random number between 0 and the lowest detected value for each metabolite

df <- subset(ft_75, SampleType != "PB")[, which(names(ft_75) %in% keep)]

has_zeros      <- apply(df, 2, function(x) any(x == 0))
zero_features  <- names(df)[has_zeros]
s_ft_has_zeros <- df[, zero_features]

lods     <- apply(s_ft_has_zeros, 2, function(x) min(x[x > 0]))
rand_imp <- sapply(1:length(lods), function(n) { runif(1000, 0, lods[n]) })

for (i in 1:length(lods)) {
  s_ft_has_zeros[, i][s_ft_has_zeros[, i] == 0] <- sample(rand_imp[, i], colSums(s_ft_has_zeros == 0)[i])
}

ft_75_imp <- cbind.data.frame(s_ft_has_zeros, df[, !colnames(df) %in% colnames(s_ft_has_zeros)])

# ---------------------------------------------------------------------------
# PCoA pre-normalisation
# ---------------------------------------------------------------------------
# We can estimate the batch effect prior batch correction using principal coordinate analysis and subsequent PERMANOVA. 
# This can be compared to the batch effect after batch correction has been applied in order to estimate the effect

# Remove Samples without EDTA or with low volume

ft_75_imp_remout <- ft_75_imp[paste("CHD",
                                     sapply(1:nrow(ft_75_imp), function(n) { strsplit(rownames(ft_75_imp), "_")[[n]][2] }),
                                     sapply(1:nrow(ft_75_imp), function(n) { strsplit(rownames(ft_75_imp), "_")[[n]][4] }), sep = "_") %in%
                                subset(md, anticoagulant_outlier == 0 & low_vol == 0)$SampleName, ]

distm <- dist(scale(ft_75_imp_remout[grep("_S", rownames(ft_75_imp_remout)), ]), method = 'euclidean')

PcoA         <- cmdscale(distm, k = 2, eig = T, add = T)
PcoA_points  <- as.data.frame(PcoA$points)
variance     <- round(PcoA$eig * 100 / sum(PcoA$eig), 1)
names(PcoA_points)[1:2] <- c('PCoA1', 'PCoA2')

identical(rownames(PcoA_points), rownames(subset(subset(ft_75[rownames(ft_75) %in% rownames(ft_75_imp_remout), ], SampleType == "S"))))

# PERMANOVA
adonres <- adonis2(distm ~ subset(ft_75[rownames(ft_75) %in% rownames(ft_75_imp_remout), ], SampleType == "S")$Batch)
adonres
p_value <- adonres$`Pr(>F)`[1]
r2      <- adonres$R2[1]

p_pcoa_batch <- ggplot(PcoA_points, aes(x = PCoA1, y = PCoA2,
                                         colour = subset(ft_75[rownames(ft_75) %in% rownames(ft_75_imp_remout), ], SampleType == "S")$Batch,
                                         label = row.names(PcoA))) +
  geom_point(size = 2.5) +
  scale_colour_jama() +
  xlab(paste('PCoA1', variance[1], '%', sep = ' ')) +
  ylab(paste('PCoA2', variance[2], '%', sep = ' ')) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  annotate("text", x = min(PcoA_points$PCoA1), y = max(PcoA_points$PCoA2),
           label = paste("p-value =", round(p_value, 4), "\nR² =", round(r2, 4)),
           hjust = 0, vjust = 1, size = 4, fontface = "bold")

p_pcoa_batch
ggsave(file.path(output_dir, "QC/pcoa_prenorm.pdf"), p_pcoa_batch)

# ---------------------------------------------------------------------------
# WaveICA normalisation
# ---------------------------------------------------------------------------
# In order to decrease the batch effects that was observed in the PERMANOVA, we utilize the correction algorithm, 
# WaveICA, designed to deal with batch effects in untargeted metabolomcis data. WaveICA can handle both intra- 
# and interbatch drift. More information about the algorithm can be found in the published https://doi.org/10.1007/s11306-021-01839-7
# and the code is also available https://github.com/dengkuistat/WaveICA_2.0

Injection_order <- subset(ft_75, SampleType != "PB")$RunOrder +
  80 * (as.numeric(as.character(subset(ft_75, SampleType != "PB")$Batch)) - 1)

ft_norm <- WaveICA2.0::WaveICA_2.0(data = ft_75_imp, Injection_Order = Injection_order, alpha = 0, Cutoff = 0.1, K = 10)

# ---------------------------------------------------------------------------
# Scale data and cap outliers
# ---------------------------------------------------------------------------
# Each metabolite is centered and univariance scaled, i.e. mean=0 and sd=1. Outliers are defined as observations 
# that deviate more than 5 sd units from the mean in either direction. The outer boundary for each metabolite 
# is set to [-5,5], so that high outliers will be converted to 5 and low outliers to -5.

ft_75_out <- scale(ft_norm$data_wave)

ft_75_out[ft_75_out > 5]  <- 5
ft_75_out[ft_75_out < -5] <- -5

ft_75_out <- cbind.data.frame(rownames(ft_75_out), subset(ft_75, SampleType != "PB")[1:3], ft_75_out)
colnames(ft_75_out)[1] <- "SampleName"

ft_75_out$SampleName <- paste("CHD",
                               sapply(1:nrow(ft_75_out), function(n) { strsplit(ft_75_out$SampleName, "_")[[n]][2] }),
                               sapply(1:nrow(ft_75_out), function(n) { strsplit(ft_75_out$SampleName, "_")[[n]][4] }), sep = "_")

# EDTA outliers and samples with recorded low volume at the laboratory work are excluded
ft_75_out[ft_75_out$SampleName %in% subset(md, anticoagulant_outlier == 0 & low_vol == 0)$SampleName,
           5:ncol(ft_75_out)] %>% head()

# ---------------------------------------------------------------------------
# PCoA post-normalisation
# ---------------------------------------------------------------------------
# Running a PCoA after the batch correction will indicate whehter WaveICA managed to 
# decrease the metabolomic differences between the batches.

distm <- dist(scale(subset(ft_75_out[ft_75_out$SampleName %in% subset(md, anticoagulant_outlier == 0 & low_vol == 0)$SampleName, ],
                            SampleType == "S")[, 5:ncol(ft_75_out)]), method = 'euclidean')

PcoA        <- cmdscale(distm, k = 2, eig = T, add = T)
PcoA_points <- as.data.frame(PcoA$points)
variance    <- round(PcoA$eig * 100 / sum(PcoA$eig), 1)
names(PcoA_points)[1:2] <- c('PCoA1', 'PCoA2')

identical(rownames(PcoA_points), rownames(subset(ft_75_out[ft_75_out$SampleName %in% subset(md, anticoagulant_outlier == 0 & low_vol == 0)$SampleName, ],
                                                  SampleType == "S")))

# PERMANOVA
adonres <- adonis2(distm ~ subset(ft_75_out[ft_75_out$SampleName %in% subset(md, anticoagulant_outlier == 0 & low_vol == 0)$SampleName, ],
                                   SampleType == "S")$Batch)
adonres
p_value <- adonres$`Pr(>F)`[1]
r2      <- adonres$R2[1]

p_pcoa_batch_postwave <- ggplot(PcoA_points, aes(x = PCoA1, y = PCoA2,
                                                   colour = subset(ft_75_out[ft_75_out$SampleName %in% subset(md, anticoagulant_outlier == 0 & low_vol == 0)$SampleName, ],
                                                                   SampleType == "S")$Batch,
                                                   label = row.names(PcoA))) +
  geom_point(size = 2.5) +
  scale_colour_jama() +
  xlab(paste('PCoA1', variance[1], '%', sep = ' ')) +
  ylab(paste('PCoA2', variance[2], '%', sep = ' ')) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  annotate("text", x = min(PcoA_points$PCoA1), y = max(PcoA_points$PCoA2),
           label = paste("p-value =", round(p_value, 4), "\nR² =", round(r2, 4)),
           hjust = 0, vjust = 1, size = 4, fontface = "bold")

p_pcoa_batch_postwave
ggsave(file.path(output_dir, "QC/pcoa_postnorm.pdf"), p_pcoa_batch_postwave)

# ---------------------------------------------------------------------------
# Batch normalisation comparison (sanity check)
# ---------------------------------------------------------------------------
# This step is a sanity check. It will compare the dataset generated by WaveICA, 
# with one that use a more simple approach to correct for batch effects, 
# by z-scoring per batch. The expected outcome is a high correlation between the two datasets

S  <- ft_75_imp
md_sub <- subset(ft_75, SampleType != "PB")[, 1:3]

batchfts   <- split(as.data.frame(S), f = md_sub$Batch)
batchft_TT <- lapply(batchfts, scale, center = T, scale = T)
batchft_TT <- do.call("rbind", batchft_TT)

# Remove outliers
batchft_TT[batchft_TT > 5]  <- 5
batchft_TT[batchft_TT < -5] <- -5
batchft_TT <- scale(batchft_TT)

ft_75_bnorm <- cbind.data.frame(md_sub, batchft_TT)
ft_75_bnorm <- cbind.data.frame(rownames(ft_75_bnorm), ft_75_bnorm)
colnames(ft_75_bnorm)[1] <- "SampleName"

dim(ft_75_bnorm)
dim(ft_75_out)

#Calculate the correlation coefficients between each metabolite in the WaveICA and batch norm datasets
suppressWarnings({
  cor_out <- sapply(5:ncol(ft_75_out), function(n) { cor(ft_75_bnorm[, n], ft_75_out[, n]) })
})

hist(cor_out)

#Expected: the correlation is high, but some metabolites are transformed a bit more with WaveICA than others.

# ---------------------------------------------------------------------------
# Binarise sparse features and remove batch-associated artifacts
# ---------------------------------------------------------------------------
# Sparse features will be converted to binary variables. These variables are less suitable for imputation but may still contain useful information.
S_25 <- subset(ft_25, SampleType == "S")
S_25[, 4:ncol(S_25)][S_25[, 4:ncol(S_25)] > 0] <- 1

ft_25_out <- cbind.data.frame(rownames(S_25),
                                S_25[, p.adjust(sapply(4:ncol(S_25), function(n) {
                                  chisq.test(S_25[, n], S_25$Batch, simulate.p.value = TRUE)$p.value
                                }), method = "fdr") > 0.05])

colnames(ft_25_out)[1] <- "SampleName"

ft_25_out$SampleName <- paste("CHD",
                               sapply(1:nrow(ft_25_out), function(n) { strsplit(rownames(ft_25_out), "_")[[n]][2] }),
                               sapply(1:nrow(ft_25_out), function(n) { strsplit(rownames(ft_25_out), "_")[[n]][4] }), sep = "_")

head(ft_25_out)

# ---------------------------------------------------------------------------
# Boxplots after correction
# ---------------------------------------------------------------------------
# Boxplots will be created for each continuous metabolite in the dataset. 
# We want the batch effect to be small (or none) and the overall variance in the analytical samples 
# to be larger than in the EC and Pools. The average in the Pools should be similar to the average 
# of the analytical samples. Could be useful to determine the analytical quality of interesting 
# biomarkers from the downstream analysis. 

pdf(file.path(output_dir, "QC/Boxplots/Boxplots_samp_postcorrection.pdf"))

plot_list <- list()
for (i in 5:ncol(ft_75_out)) {
  plot_list[[i]] <- ggplot(ft_75_out, aes(x = SampleType, y = ft_75_out[, i], fill = Batch)) +
    geom_boxplot() +
    theme_minimal() +
    scale_fill_jama() +
    ylab("") +
    ggtitle(paste0("feature", names(ft_75_out)[i]))
  print(plot_list[[i]])
}

dev.off()

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

write.csv(ft_75_out, file.path(output_dir, "chd_quant_all_cont_blankrm_wavenorm.csv"))

sessionInfo()