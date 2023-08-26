# # A script to read in quasi alignments with Salmon and use EdgeR to assess differential gene expression
# BiocManager::install("tximport")
# BiocManager::install("edgeR")
# BiocManager::install("BUSpaRse")

library(tidyverse) # Tidyverse for handling large datasets
library(tximport) # Tximport for reading the salmon output files
library(edgeR) # EdgeR for differential expression
library(pheatmap) # Pretty heatmaps
library(ggbeeswarm) # This is for plotting nicer scatterplots across groups
library(fishualize) # Get nice fish-based colours
library(rjson) # Read json files to get mapping percentages from salmon

# Read in the gene to transcript map. Switch the order of columns to match the order tximport wants
gene_trans_map <- read_tsv("../Trinity_gill_genetransmap.txt", col_names = c("gene_id", "transcript_id")) %>% 
  select(transcript_id, gene_id)

# Read in the metadata, filter for the samples which are present (not all were sequenced from the original experiment)
metadata <- readxl::read_excel("../Luke's Gill samples_Jan3_2022.xlsx", sheet = "metadata_reformat")
metadata_filtered <- metadata %>% 
  filter(!is.na(Data_ID)) %>% 
  mutate(Treatment_Time = paste0(Treatment, "_", Hours)) # Add a column with all the treatment information in it, to simplify contrasts for a later EdgeR model

# Follow the tximport vignette for reading in salmon quant.sf files 
# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon
# Get a list of the file paths for each sample
files <- file.path("salmon_quant", metadata_filtered$Data_ID, "quant.sf")
names(files) <- paste0(c(metadata_filtered$Data_ID))

# Import the salmon quant.sf files. Run this line for 'gene'-level counts using the gene transcript map
#txi.salmon <- tximport(files, type = "salmon", tx2gene = gene_trans_map)

# Run this line for 'transcript' level counts that do not summarize data to the 'gene' level
txi.salmon <- tximport(files, type = "salmon", txOut = TRUE)

head(txi.salmon$counts)

# Convert the tximport-formatted data to the edgeR format following the tximport vignette
# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#edgeR
counts <- txi.salmon$counts
normMat <- txi.salmon$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- counts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
edgeR_data <- DGEList(counts)

# Input the offset matrix calculated previously
edgeR_data <- scaleOffset(edgeR_data, normMat)

# filtering using the design information
#design <- model.matrix(~ Hours * Temp * PCO2, data = metadata_filtered)
#design <- model.matrix(~ Treatment, data = metadata_filtered)
# Rename columns to keep R happy, make contrasts for the interaction model
#colnames(design) <- make.names(c("Intercept", "Hours", "Temp", "PCO2", "Hours:Temp","Hours:PCO2","Temp:PCO2","Hours:Temp:PCO2"))
#interaction_contrast <- makeContrasts(hours = Hours - Temp, levels = design)

# A design matrix based on time series. Following section 4.8 "Time course RNA-seq experiments of Drosophila melanogaster" from the EdgeR user's guide
# https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#poly_hours <- poly(metadata_filtered$Hours, degree = 2)
#design_time <- model.matrix(~ poly_hours)

#keep <- filterByExpr(edgeR_data, design_time)

# Filter out genes without expression in any of the groups
#edgeR_data <- edgeR_data[keep, ]
# edgeR_data is now ready for estimate dispersion functions see edgeR User's Guide
#nrow(edgeR_data$counts)

# Estimate dispersion
#edgeR_data <- estimateDisp(edgeR_data, design_time)
#plotBCV(edgeR_data)
#sqrt(edgeR_data$common.dispersion)

# Take a look at a multidimensional scaling plot
#plotMDS(edgeR_data, top = nrow(edgeR_data$counts), labels = c(paste0(metadata_filtered$Treatment, "_", metadata_filtered$Hours)))

# Fit the GLM using the edgeR data and the design
#edger_fit <- glmQLFit(edgeR_data, design_time, robust = TRUE)
#plotQLDisp(edger_fit)

# Conduct an F-test on 2 degrees of freedom for each gene
#edger_fit <- glmQLFTest(edger_fit, coef=2:3)

# What are the 30 most significant genes?
#top_tags <- as.data.frame(topTags(edger_fit, n=30))
#top_tags

# How many genes were significant and not significant at q < 0.05?
#summary(decideTests(edger_fit))

# We visualize the fitted spline curves for the top four genes. We start by computing the observed and fitted log-CPM values for each gene
#logCPM.obs <- cpm(edgeR_data, log = TRUE, prior.count = edger_fit$prior.count)
#logCPM.fit <- cpm(edger_fit, log=TRUE)

# Plot the log2 CPM values for the genes that changed most by time
# par(mfrow=c(2,2))
# for(i in 1:4) {
#   gene_id <- row.names(top_tags)[i]
#   logCPM.obs.i <- logCPM.obs[gene_id,]
#   logCPM.fit.i <- logCPM.fit[gene_id,]
#   plot(metadata_filtered$Hours, logCPM.obs.i, ylab="log-CPM", main = gene_id, pch=16)
#   lines(metadata_filtered$Hours, logCPM.fit.i, col="red", lwd=2)
#   }
# dev.off() # Just to reset par

# Re-doing the EdgeR analysis with a different design formula, using a combined measure of time, temperature, and PCO2 to simplify interpretation of contrasts
# Use a 0 in the design formula to not to include an intercept column and instead to include a column for each group
design <- model.matrix(~ 0 + Treatment_Time, data = metadata_filtered) 

# Turn the design matrix into a tibble to remove the redundant 'treatment_time' string from each column name
design <- as_tibble(design) %>% 
  rename_with(~str_remove(., 'Treatment_Time'))

# Turn the design tibble back into a matrix
design <- as.matrix(design)

# Use make.names to take out characters that R does not like in column names
colnames(design) <- make.names(c(colnames(design)))

# Creating a DGEList object for use in edgeR.
edgeR_data_treatment <- DGEList(counts)

# Input the offset matrix calculated previously
edgeR_data_treatment <- scaleOffset(edgeR_data_treatment, normMat)

# Filter genes by abundance
keep_treatment <- filterByExpr(edgeR_data_treatment, design)
# How many are we keeping?
sum(keep_treatment) # 162,557

# Filter out genes without expression in any of the groups
edgeR_data_treatment <- edgeR_data_treatment[keep_treatment, ]
# edgeR_data is now ready for estimate dispersion functions see edgeR User's Guide

# Estimate dispersion
edgeR_data_treatment <- estimateDisp(edgeR_data_treatment, design)
plotBCV(edgeR_data_treatment)
sqrt(edgeR_data_treatment$common.dispersion)

# Take a look at a multidimensional scaling plot
plotMDS(edgeR_data_treatment, top = nrow(edgeR_data_treatment$counts), labels = c(metadata_filtered$Treatment_Time))

# Create a dataframe from the MDS to plot with ggplot. Note that  because all genes are used in the plotMDS function, it is *not* an MDS but is instead a PCA. See the plotMDS documentation for details, but all transcripts were used here with top = nrow(), thus it is not an MDS which would have used some subset of transcripts
#MDS_data <- plotMDS(edgeR_data_treatment, top = nrow(edgeR_data_treatment$counts), labels = c(metadata_filtered$Treatment_Time), plot = FALSE)

# Here, an MDS is run with the top 10% of transcripts (16,256) by leading log2-fold changes
MDS_data <- plotMDS(edgeR_data_treatment, top = nrow(edgeR_data_treatment$counts)*0.1, labels = c(metadata_filtered$Treatment_Time), plot = FALSE)
MDS_data$var.explained # Look at variance explained. 5.67% for dimension 1, 3.23% for dimension 2
MDS_data_tibble <- tibble(x_coord = MDS_data$x, y_coord = MDS_data$y, Time = as.factor(metadata_filtered$Hours), Treatment = metadata_filtered$Treatment)

plotMDS(edgeR_data_treatment, top = nrow(edgeR_data_treatment$counts)*0.1, labels = c(metadata_filtered$Treatment_Time), plot = TRUE)
# Plot the MDS
MDS_plot <- ggplot(data = MDS_data_tibble, aes(x = x_coord, y = y_coord, color = Treatment, shape = Time)) +
  geom_point(size = 4.5, alpha = 0.9) +
  #xlab(bquote(Leading~log[2]*"-fold change Dimension 1 (5.6%)")) + # The PCA is not an MDS, so the leading log2 fold change line is wrong
  xlab(bquote("Leading 10% log"[2]*"-Fold Change Dimension 1 (5.67%)")) +
  ylab(bquote("Leading 10% log"[2]*"-Fold Change Dimension 2 (3.23%)")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.85, 0.78), legend.background = element_blank(), legend.key = element_blank(), legend.key.size = unit(0.4, 'cm'), legend.spacing.y = unit(0.1, 'cm'), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
ggsave(filename = "edgeR_MDS.pdf", plot = MDS_plot, dpi = 3000, scale = 1.34)

# Fit the GLM using the edgeR data and the design
edger_fit_treatment <- glmQLFit(edgeR_data_treatment, design, robust = TRUE)
plotQLDisp(edger_fit_treatment)

# Define contrasts for the design with treatment and time.
# First, design a series of control contrasts. These will be used to subtract out genes that changed in the control fish over time from the PCO2 and Temp comparisons.
treatment_contrast <- makeContrasts(control_6v0 = Control_6 - Control_0, 
                                    control_168v6 = Control_168 - Control_6,
                                    control_168v0 = Control_168 - Control_0,
                                    pCO2_6v0 = PCO2_6 - PCO2_0,
                                    pCO2_168v6 = PCO2_168 - PCO2_6,
                                    pCO2_168v0 = PCO2_168 - PCO2_0,
                                    temp_6v0 = Temp_6 - Temp_0,
                                    temp_168v6 = Temp_168 - Temp_6,
                                    temp_168v0 = Temp_168 - Temp_0,
                                    pCO2.temp_6v0 = Temp.PCO2_6 - Temp.PCO2_0,
                                    pCO2.temp_168v6 = Temp.PCO2_168 - Temp.PCO2_6,
                                    pCO2.temp_168v0 = Temp.PCO2_168 - Temp.PCO2_0,
                                    levels = design)
treatment_contrast

# Pulling results for the control contrasts
quasi_test_control6v0 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"control_6v0"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_control6v0, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_control6v0)

quasi_test_control168v6 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"control_168v6"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_control168v6, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_control168v6)

quasi_test_control168v0 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"control_168v0"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_control168v0, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_control168v0)

# Pull a table of genes that were significant in any of the control treatments. Put them all into one big table with a new column for which contrast the gene is from. Make sure to keep gene information. 
# 6,747 unique genes were significant in these control comparisons
sig_control_genes <- as_tibble(bind_rows(as.data.frame(topTags(quasi_test_control6v0, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
      mutate(contrast = "control_6v0") %>% 
      rownames_to_column(var = "gene"), 
      as.data.frame(topTags(quasi_test_control168v6, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
      mutate(contrast = "control_168v6") %>% 
      rownames_to_column(var = "gene"), 
      as.data.frame(topTags(quasi_test_control168v0, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
      mutate(contrast = "control_168v0") %>% 
      rownames_to_column(var = "gene")))
nrow(distinct(sig_control_genes, gene))


# Pull results for the within-PCO2 contrasts
quasi_test_pco2.6v0 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"pCO2_6v0"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_pco2.6v0, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_pco2.6v0)

quasi_test_pco2.168v6 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"pCO2_168v6"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_pco2.168v6, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_pco2.168v6)

quasi_test_pco2.168v0 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"pCO2_168v0"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_pco2.168v0, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_pco2.168v0)

# Pull a table of genes significant in any of the within-pCO2 treatments. Collate them into a table with contrast information, and keep gene names.
# 16,399 unique genes are significant in these pCO2 contrasts.
sig_pco2_genes <- as_tibble(bind_rows(as.data.frame(topTags(quasi_test_pco2.6v0, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
                      mutate(contrast = "pCO2_6v0") %>% 
                      rownames_to_column(var = "gene"), 
                    as.data.frame(topTags(quasi_test_pco2.168v6, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
                      mutate(contrast = "pCO2_168v6") %>% 
                      rownames_to_column(var = "gene"), 
                    as.data.frame(topTags(quasi_test_pco2.168v0, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
                      mutate(contrast = "pCO2_168v0") %>% 
                      rownames_to_column(var = "gene")))
nrow(distinct(sig_pco2_genes, gene))

# Pull results for the within-temperature contrasts
quasi_test_temp.6v0 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"temp_6v0"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_temp.6v0, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_temp.6v0)

quasi_test_temp.168v6 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"temp_168v6"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_temp.168v6, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_temp.168v6)

quasi_test_temp.168v0 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"temp_168v0"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_temp.168v0, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_temp.168v0)

# Pull a table of genes significant in any of the within-temperature treatments. Collate them into a table with contrast information, and keep gene names.
# 17,480 unique genes are significant in these temperature contrasts.
sig_temp_genes <- as_tibble(bind_rows(as.data.frame(topTags(quasi_test_temp.6v0, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
                                        mutate(contrast = "temp_6v0") %>% 
                                        rownames_to_column(var = "gene"), 
                                      as.data.frame(topTags(quasi_test_temp.168v6, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
                                        mutate(contrast = "temp_168v6") %>% 
                                        rownames_to_column(var = "gene"), 
                                      as.data.frame(topTags(quasi_test_temp.168v0, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
                                        mutate(contrast = "temp_168v0") %>% 
                                        rownames_to_column(var = "gene")))
nrow(distinct(sig_temp_genes, gene))

# Pull results for the within-pCO2 and temperature contrasts
quasi_test_pco2.temp.6v0 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"pCO2.temp_6v0"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_pco2.temp.6v0, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_pco2.temp.6v0)

quasi_test_pco2.temp.168v6 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"pCO2.temp_168v6"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_pco2.temp.168v6, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_pco2.temp.168v6)

quasi_test_pco2.temp_168v0 <- glmQLFTest(edger_fit_treatment, contrast = treatment_contrast[,"pCO2.temp_168v0"])
# Look at how many genes were significant at q < 0.05, return as many genes as are significant with n
nrow(topTags(quasi_test_pco2.temp_168v0, p.value = 0.05, n = nrow(edger_fit_treatment)))
# Plot the results for the contrast
plotMD(quasi_test_pco2.temp_168v0)

# Pull a table of genes significant in any of the within-pCO2 and temperature treatments. Collate them into a table with contrast information, and keep gene names.
# 17,494 unique genes are significant in these pCO2 and temperature contrasts.
sig_pCO2_temp_genes <- as_tibble(bind_rows(as.data.frame(topTags(quasi_test_pco2.temp.6v0, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
                                        mutate(contrast = "pCO2_temp_6v0") %>% 
                                        rownames_to_column(var = "gene"), 
                                      as.data.frame(topTags(quasi_test_pco2.temp.168v6, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
                                        mutate(contrast = "pCO2_temp_168v6") %>% 
                                        rownames_to_column(var = "gene"), 
                                      as.data.frame(topTags(quasi_test_pco2.temp_168v0, p.value = 0.05, n = nrow(edger_fit_treatment))) %>% 
                                        mutate(contrast = "pCO2_temp_168v0") %>% 
                                        rownames_to_column(var = "gene")))
nrow(distinct(sig_pCO2_temp_genes, gene))

# Filter out the genes significant within the control contrasts. Do this by taking the intersection of the control and contrast genes to find matching genes, then filtering them from the results tables.
sig_pco2_genes <- filter(sig_pco2_genes, !(gene %in% c(intersect(sig_control_genes$gene, sig_pco2_genes$gene))))
nrow(distinct(sig_pco2_genes, gene))

sig_temp_genes <- filter(sig_temp_genes, !(gene %in% c(intersect(sig_control_genes$gene, sig_temp_genes$gene))))
nrow(distinct(sig_temp_genes, gene))

sig_pCO2_temp_genes <- filter(sig_pCO2_temp_genes, !(gene %in% c(intersect(sig_control_genes$gene, sig_pCO2_temp_genes$gene))))
nrow(distinct(sig_pCO2_temp_genes, gene))

# After filtering out genes that changed within the control group, significant genes that remained were:
# 13,033 for pCO2
# 13,243 for temp
# 13,943 for pCO2 and temp combined

# Write out tables for DGE with the control genes subtracted, but *not* unique genes
write_tsv(x = sig_pco2_genes, file = "sig_pCO2_DGE_nonunique.txt")
write_tsv(x = sig_temp_genes, file = "sig_temp_DGE_nonunique.txt")
write_tsv(x = sig_pCO2_temp_genes, file = "sig_pCO2_temp_DGE_nonunique.txt")


# Subtract out genes from the different contrasts to get 'unique' genes in each contrast. 
# First, pCO2-unique genes gotten by subtracting out genes that appear in the temperature results
unique_pco2_genes <- filter(sig_pco2_genes, !(gene %in% c(intersect(sig_temp_genes$gene, sig_pco2_genes$gene))))
nrow(distinct(unique_pco2_genes, gene))

# Temperature-specific genes
unique_temp_genes <- filter(sig_temp_genes, !(gene %in% c(intersect(sig_pco2_genes$gene, sig_temp_genes$gene))))
nrow(distinct(unique_temp_genes, gene))

# pCO2 and temperature combined, unique genes
unique_pco2_temp_genes <- filter(sig_pCO2_temp_genes, !(gene %in% c(intersect(sig_pco2_genes$gene, sig_pCO2_temp_genes$gene))), !(gene %in% c(intersect(sig_temp_genes$gene, sig_pCO2_temp_genes$gene))))
nrow(distinct(unique_pco2_temp_genes, gene))

# After filtering out genes in the other contrasts to get 'unique' genes for each treatment. The genes that remained were:
# 9,414 for pCO2
# 9,264 for temp
# 7,007 for pCO2 and temp combined

# Write out the tables of unique genes for enrichment analyses
write_tsv(x = unique_pco2_genes, file = "unique_pCO2_genes.txt")
write_tsv(x = unique_temp_genes, file = "unique_temp_genes.txt")
write_tsv(x = unique_pco2_temp_genes, file = "unique_pCO2_temp_genes.txt")


# Dig into the 168 v 6 hour contrasts. The goal here is to find acute vs 'chronic' (or sub-chronic) responses to the different treatments
pco2_168v6 <- filter(unique_pco2_genes, contrast == "pCO2_168v6")
nrow(distinct(pco2_168v6, gene))

temp_168v6 <- filter(unique_temp_genes, contrast == "temp_168v6")
nrow(distinct(temp_168v6, gene))

pco2_temp_168v6 <- filter(unique_pco2_temp_genes, contrast == "pCO2_temp_168v6")
nrow(distinct(pco2_temp_168v6, gene))

# Within the 168 v 6 hour contrasts, the unique genes that remained were:
# 5,056 in pCO2
# 7,490 in temp
# 5,041 in pCO2 and temp combined

# Write out these tables of 168 v 6 hour results
write_tsv(x = pco2_168v6, file = "pco2_168v6.txt")
write_tsv(x = temp_168v6, file = "temp_168v6.txt")
write_tsv(x = pco2_temp_168v6, file = "pco2_temp_168v6.txt")

##################################################################################################
# Create a heatmap
# Pull a counts per million matrix for making a heatmap 
counts_per_million <- edgeR::cpm(edgeR_data, normalized.lib.sizes = TRUE, offset = edgeR_data$offset, log = TRUE)

# Convert the counts per million matrix into a tibble, replace the | character with _ to allow for selecting genes
counts_tibble <- counts_per_million %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "gene") %>% 
  mutate(gene = str_replace(gene, pattern = "\\|", replacement = "_")) %>% 
  column_to_rownames(., var = "gene")

# Select all unique genes significant within the different treatments. Replace | with _
select <- unique(c(unique_pco2_genes$gene, unique_temp_genes$gene, unique_pco2_temp_genes$gene)) %>% 
  str_replace(pattern = "\\|", replacement = "_")
length(select)

# Re-order the metadata by treatment and then time (as opposed to time and then treatment)
metadata_treatment_reorder <- metadata_filtered %>% 
  arrange(Treatment, Hours) %>% 
  mutate(Treatment_num = as.numeric(as.factor(Treatment)), ID = Data_ID, Hours = as.factor(Hours)) %>% 
  column_to_rownames(var = "Data_ID")

sample_col <- select(metadata_treatment_reorder, Hours, Treatment) %>% 
  mutate(Treatment = case_when(Treatment == "Control" ~ "Control",
                               Treatment == "PCO2" ~ "pCO2",
                               Treatment == "Temp" ~ "Temp",
                               Treatment == "Temp+PCO2" ~ "Temp+pCO2"))

# Look at 4 fishualize colours for use with the heatmap
fish(n = 4, option = "Lampris_guttatus")

# Pull counts for all these genes, make a heatmap. Write the heatmap to a pdf or png (change the file extension to change format)
pheatmap(counts_tibble[select, metadata_treatment_reorder$ID] %>% drop_na(), cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, show_colnames = FALSE, annotation_col = sample_col, annotation_colors = list(Hours = c("0" = "#fdae61", "6" = "#d73027", "168" = "#2c7bb6"), Treatment = c("Control" = "#13317DFF", "pCO2" = "#7CC2D7FF", "Temp" = "#C1AF95FF", "Temp+pCO2" = "#F75107FF")), filename = "uniquetranscripts_heatmap.png", fontsize_number = 18)



# Try a gg version of the heatmap. It ended up not being used because adding labels in ggplot did not look as good as pheatmap

#test <- pheatmap(counts_tibble[select, metadata_treatment_reorder$ID] %>% drop_na(), cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, show_colnames = FALSE, annotation_col = sample_col, annotation_colors = list(Hours = c("0" = "#fdae61", "6" = "#d73027", "168" = "#2c7bb6"), Treatment = c("Control" = "#13317DFF", "pCO2" = "#7CC2D7FF", "Temp" = "#C1AF95FF", "Temp+pCO2" = "#F75107FF")), plot = FALSE)

# heatmap_tibble <- counts_tibble[select, metadata_treatment_reorder$ID] %>% 
#   drop_na() %>% 
#   rownames_to_column(var = "transcript_ID") %>% 
#   pivot_longer(cols = L61:L88) %>% 
#   mutate(transcript_ID = factor(transcript_ID, levels = c(test$tree_row$labels)),
#          name = factor(name, levels = c(metadata_treatment_reorder$ID)))
# 
# test_plot <- ggplot(heatmap_tibble, aes(x = name, y = transcript_ID, fill = value)) + 
#   geom_tile()
# 
# ggsave(filename = "ggheatmap.pdf", plot = test_plot, width = 20)

#pheatmap(counts_tibble[c(unique(unique_pco2_genes$gene)) %>% str_replace(pattern = "\\|", replacement = "_"),] %>% drop_na(), cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, main = "pCO2 Unique Heatmap")

#pco2_select <- c(unique(unique_pco2_genes$gene))

# Write out the log2-transformed counts per million matrix for the WGCNA
write_tsv(x = as_tibble(counts_per_million, rownames = NA) %>% rownames_to_column(var = "Gene_ID"), file = "log2_cpm.txt")
##################################################################################################
# Pull CPM values for plotting individual genes of interest
counts_per_million_nolog <- edgeR::cpm(edgeR_data, normalized.lib.sizes = TRUE, offset = edgeR_data$offset, log = FALSE)

# Write out the non-logged CPM values for gene correlations
write_tsv(x = as_tibble(counts_per_million_nolog, rownames = NA) %>% rownames_to_column(var = "Gene_ID"), file = "cpm_nolog.txt")

# Pull the counts for counts of atp1a2
atp1a2_cpm <- as_tibble(counts_per_million_nolog["TR16620|c0_g1_i1",], rownames = NA) %>% 
  rownames_to_column(var = "Data_ID") %>% 
  left_join(., metadata_filtered) %>% 
  rename(ATP1A2_cpm = value)
atp1a2_cpm$Hours <- factor(atp1a2_cpm$Hours, levels = c("0", "6", "168"))


ggplot(data = atp1a2_cpm, aes(x = Hours, y = ATP1A2_cpm, group = Treatment)) + 
  geom_beeswarm(size = 2) +
  facet_wrap(~ Treatment) +
  ylab(bquote(italic("atp1a2")~Counts~per~Million)) +
  theme_bw(base_size = 18)

# Pull the counts for counts of tbx1
# Both TR75533|c3_g1 and TR47161|c0_g2 were annotated to tbx1
tbx1_cpm <- as_tibble(counts_per_million_nolog["TR47161|c0_g2_i1",], rownames = NA) %>% 
  rownames_to_column(var = "Data_ID") %>% 
  left_join(., metadata_filtered) %>% 
  rename(TBX1_cpm = value)
tbx1_cpm$Hours <- factor(tbx1_cpm$Hours, levels = c("0", "6", "168"))

ggplot(data = tbx1_cpm, aes(x = Hours, y = TBX1_cpm, group = Treatment)) + 
  geom_beeswarm(size = 2) +
  facet_wrap(~ Treatment) +
  ylab(bquote(italic("tbx1")~Counts~per~Million)) +
  theme_bw(base_size = 18)

# Look at SLC6A4 transcripts. These are the interesting ones
TR49140|c2_g1_i8
TR49140|c2_g1_i3
TR49140|c2_g1_i4 # This is the one that was actually significant in the DGE comparisons
slc6a4_cpm <- as_tibble(counts_per_million_nolog["TR49140|c2_g1_i4",], rownames = NA) %>% 
  rownames_to_column(var = "Data_ID") %>% 
  left_join(., metadata_filtered) %>% 
  rename(cpm = value)
slc6a4_cpm$Hours <- factor(slc6a4_cpm$Hours, levels = c("0", "6", "168"))

# Plot the slc6a4 CPM in the different treatments
slc6a4_plot <- ggplot(data = slc6a4_cpm, aes(x = Hours, y = cpm, group = Treatment)) + 
  geom_violin(aes(group = Hours)) +
  geom_point(size = 2, fill = "white", shape = 1) +
  facet_wrap(~ Treatment) +
  ylab(bquote(italic("slc6a4")~Counts~per~Million)) +
  theme_bw(base_size = 18)

# Save the slc6a4 plot as a supplemental figure
ggsave(filename = "slc6a4_cpm.pdf", plot = slc6a4_plot, dpi = 2000)

# Look at Na/K ATPase transcripts
TR36101|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-1-interacting protein 3
TR41650|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1
TR41650|c0_g2_i1 # Sodium/potassium-transporting ATPase subunit alpha-1
TR41650|c0_g3_i1 # Sodium/potassium-transporting ATPase subunit alpha-1
TR22292|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1, seemed to respond to pCO2 only?
TR22292|c0_g1_i3 # Sodium/potassium-transporting ATPase subunit alpha-1
TR22292|c1_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1
TR22292|c1_g1_i2 # Sodium/potassium-transporting ATPase subunit alpha-1
TR22292|c1_g1_i3 # Sodium/potassium-transporting ATPase subunit alpha-1
TR22292|c1_g1_i4 # Sodium/potassium-transporting ATPase subunit alpha-1
TR22292|c1_g1_i5 # Sodium/potassium-transporting ATPase subunit alpha-1
TR22292|c1_g1_i7 # Sodium/potassium-transporting ATPase subunit alpha-3
TR45365|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-2, decreasing with control and temp, increasing with combined temp+pCO2?
TR45365|c0_g1_i2 # Sodium/potassium-transporting ATPase subunit beta-2, increasing with pCO2 only, decreasing with control?
TR45365|c0_g1_i3 # Sodium/potassium-transporting ATPase subunit beta-2
TR45365|c0_g1_i4 # Sodium/potassium-transporting ATPase subunit beta-2
TR16620|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-2, decreasing with pCO2 only, slight increase in temp and temp+pCO2?
TR55914|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha
TR85165|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1, increasing with control and temp+pCO2?
TR74551|c2_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-1, increasing with control and temp+pCO2, slight increase with temp?
TR74551|c2_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-1, increasing with control and temp+pCO2, slight increase with temp?
TR74551|c2_g1_i2 # Sodium/potassium-transporting ATPase subunit beta-1, increasing with temp and temp+pCO2? MIGHT MATCH ENZYME ACTIVITY
TR74551|c2_g1_i3 # Sodium/potassium-transporting ATPase subunit beta-1, increasing with control and temp+pCO2, slight increase with temp?
TR74551|c2_g1_i4 # Sodium/potassium-transporting ATPase subunit beta-1, increasing with temp and temp+pCO2? MIGHT MATCH ENZYME ACTIVITY
TR74551|c2_g1_i5 # Sodium/potassium-transporting ATPase subunit beta-1, increasing with temp and temp+pCO2? MIGHT MATCH ENZYME ACTIVITY
TR74551|c2_g1_i6 # Sodium/potassium-transporting ATPase subunit beta-233, increasing with temp and temp+pCO2? Less than other transcripts
TR74551|c2_g1_i9 # Sodium/potassium-transporting ATPase subunit beta-1, increasing with temp and temp+pCO2? Less than other transcripts
TR36101|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-1-interacting protein 3
TR40064|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1, increase with temp, slight decrease in pCO2?
TR40064|c1_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-3
TR40580|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1, decrease with pCO2, slight decrease with temp?
TR76055|c2_g1_i2 # Sodium/potassium-transporting ATPase subunit alpha-3, slight decrease in control, increase in pCO2 and in temp at 7 days?
TR76055|c2_g1_i3 # Sodium/potassium-transporting ATPase subunit alpha-3
TR70983|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-3
TR70983|c0_g1_i2 # Sodium/potassium-transporting ATPase subunit beta-3
TR70983|c0_g1_i3 # Sodium/potassium-transporting ATPase subunit beta-3
TR76055|c2_g1_i2 # Sodium/potassium-transporting ATPase subunit alpha-3
TR76055|c2_g1_i3 # Sodium/potassium-transporting ATPase subunit alpha-3
TR70983|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-3
TR70983|c0_g1_i3 # Sodium/potassium-transporting ATPase subunit beta-3
TR41650|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1
TR41650|c0_g2_i1 # Sodium/potassium-transporting ATPase subunit alpha-1
TR41650|c0_g3_i1 # Sodium/potassium-transporting ATPase subunit alpha-1
TR16620|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-2
TR36101|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-1-interacting protein 3
TR45365|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-2, decrease in control, increase in temp+pCO2?
TR45365|c0_g1_i2 # Sodium/potassium-transporting ATPase subunit beta-2, increase in pCO2 only?
TR45365|c0_g1_i3 # Sodium/potassium-transporting ATPase subunit beta-2
TR45365|c0_g1_i4 # Sodium/potassium-transporting ATPase subunit beta-2, increase in control and temp+pCO2, decrease in pCO2 and in temp?
TR85165|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1, increase in temp+pCO2?
TR22292|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1, decrease in pCO2?
TR22292|c0_g1_i3 # Sodium/potassium-transporting ATPase subunit alpha-1, decrease in pCO2, control, and temp?
TR22292|c1_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1, increase in temp+pCO2?
TR22292|c1_g1_i2 # Sodium/potassium-transporting ATPase subunit alpha-1, increase in temp+pCO2?
TR22292|c1_g1_i3 # Sodium/potassium-transporting ATPase subunit alpha-1, decrease in pCO2?
TR22292|c1_g1_i4 # Sodium/potassium-transporting ATPase subunit alpha-1, slight increase in temp+pCO2?
TR22292|c1_g1_i5 # Sodium/potassium-transporting ATPase subunit alpha-1
TR22292|c1_g1_i7 # Sodium/potassium-transporting ATPase subunit alpha-3
TR55914|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha
TR74551|c2_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-1, increase in control and temp+pCO2?
TR74551|c2_g1_i1 # Sodium/potassium-transporting ATPase subunit beta-1, increase in control and temp+pCO2?
TR74551|c2_g1_i2 # Sodium/potassium-transporting ATPase subunit beta-1, large increase in temp+pCO2?
TR74551|c2_g1_i3 # Sodium/potassium-transporting ATPase subunit beta-1, increase in temp+pCO2?
TR74551|c2_g1_i3 # Sodium/potassium-transporting ATPase subunit beta-1, increase in control and temp+pCO2?
TR74551|c2_g1_i4 # Sodium/potassium-transporting ATPase subunit beta-1, increase in temp+pCO2, slight increase in control?
TR74551|c2_g1_i5 # Sodium/potassium-transporting ATPase subunit beta-1, increase in temp+pCO2, slight increase in control?
TR74551|c2_g1_i6 # Sodium/potassium-transporting ATPase subunit beta-233, increase in temp, slight increase in control and temp+pCO2?
TR74551|c2_g1_i9 # Sodium/potassium-transporting ATPase subunit beta-1, slight increase in temp?
TR40064|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1, increase in temp, decrease in pCO2?
TR40064|c1_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-3, slight increase in temp and temp+pCO2? MATCH TO ENZYME ACTIVITY?
TR40580|c0_g1_i1 # Sodium/potassium-transporting ATPase subunit alpha-1, decrease in control and pCO2?

as_tibble(counts_per_million_nolog["TR74551|c2_g1_i5",], rownames = NA) %>% 
  rownames_to_column(var = "Data_ID") %>% 
  left_join(., metadata_filtered) %>% 
  rename(cpm = value) %>% 
  mutate(Hours = factor(Hours, levels = c("0", "6", "168"))) %>% 
ggplot(data = ., aes(x = Hours, y = cpm, group = Treatment)) + 
  geom_violin(aes(group = Hours)) +
  geom_point(size = 2, fill = "white", shape = 1) +
  facet_wrap(~ Treatment) +
  ylab(bquote(italic("NKA")~Counts~per~Million)) +
  theme_bw(base_size = 18)

# Create a vector of NKA transripts
nka_transcripts <- c("TR55914|c0_g1_i1",
                     "TR55914|c0_g1_i1",
                     "TR41650|c0_g1_i1",
                     "TR41650|c0_g2_i1",
                     "TR41650|c0_g3_i1",
                     "TR22292|c0_g1_i1",
                     "TR22292|c0_g1_i3",
                     "TR22292|c1_g1_i1",
                     "TR22292|c1_g1_i2",
                     "TR22292|c1_g1_i3",
                     "TR22292|c1_g1_i4",
                     "TR22292|c1_g1_i5",
                     "TR85165|c0_g1_i1",
                     "TR40064|c0_g1_i1",
                     "TR40580|c0_g1_i1",
                     "TR41650|c0_g1_i1",
                     "TR41650|c0_g2_i1",
                     "TR41650|c0_g3_i1",
                     "TR85165|c0_g1_i1",
                     "TR22292|c0_g1_i1",
                     "TR22292|c0_g1_i3",
                     "TR22292|c1_g1_i1",
                     "TR22292|c1_g1_i2",
                     "TR22292|c1_g1_i3",
                     "TR22292|c1_g1_i4",
                     "TR22292|c1_g1_i5",
                     "TR40064|c0_g1_i1",
                     "TR40580|c0_g1_i1",
                     "TR16620|c0_g1_i1",
                     "TR16620|c0_g1_i1",
                     "TR22292|c1_g1_i7",
                     "TR40064|c1_g1_i1",
                     "TR76055|c2_g1_i2",
                     "TR76055|c2_g1_i3",
                     "TR76055|c2_g1_i2",
                     "TR76055|c2_g1_i3",
                     "TR22292|c1_g1_i7",
                     "TR40064|c1_g1_i1",
                     "TR74551|c2_g1_i1",
                     "TR74551|c2_g1_i1",
                     "TR74551|c2_g1_i2",
                     "TR74551|c2_g1_i3",
                     "TR74551|c2_g1_i4",
                     "TR74551|c2_g1_i5",
                     "TR74551|c2_g1_i9",
                     "TR74551|c2_g1_i1",
                     "TR74551|c2_g1_i1",
                     "TR74551|c2_g1_i2",
                     "TR74551|c2_g1_i3",
                     "TR74551|c2_g1_i3",
                     "TR74551|c2_g1_i4",
                     "TR74551|c2_g1_i5",
                     "TR74551|c2_g1_i9",
                     "TR45365|c0_g1_i1",
                     "TR45365|c0_g1_i2",
                     "TR45365|c0_g1_i3",
                     "TR45365|c0_g1_i4",
                     "TR45365|c0_g1_i1",
                     "TR45365|c0_g1_i2",
                     "TR45365|c0_g1_i3",
                     "TR45365|c0_g1_i4",
                     "TR74551|c2_g1_i6",
                     "TR74551|c2_g1_i6",
                     "TR70983|c0_g1_i1",
                     "TR70983|c0_g1_i2",
                     "TR70983|c0_g1_i3",
                     "TR70983|c0_g1_i1",
                     "TR70983|c0_g1_i3")

# Pull counts for the transcripts
nka_counts <- as_tibble(counts_per_million_nolog[nka_transcripts,], rownames = NA) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Data_ID") %>% 
  left_join(., metadata_filtered, by = "Data_ID") %>% 
  select(-`...8`, -`...9`, )
  #rename(gene = transcript) %>% 
  mutate(Hours = factor(Hours, levels = c("0", "6", "168")))
############################################################################################
# Create an empty tibble with two columns- individual ID and an empty column for mapping percent
mapping_result <- tibble(ID = metadata_filtered$Data_ID, mapping_perc = NA) %>% 
  arrange(ID)

# Loop through the Salmon outputs by individual ID and pull mapping results from the meta info files
for (sample in 1:nrow(mapping_result)){
  json_file <- fromJSON(paste(readLines(paste0("salmon_quant/", mapping_result$ID[sample], "/aux_info/meta_info.json")), collapse=""))
  mapping_result$mapping_perc[sample] <- json_file$percent_mapped
}

# Write out this file with mapping results to attach to the other sample data
write_tsv(x = mapping_result, file = "mapping_percents.txt")
