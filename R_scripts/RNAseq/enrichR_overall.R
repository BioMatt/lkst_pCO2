library(tidyverse)
library(enrichR)

# A function to take an EnrichrR table and split the GO term and definition for input into Revigo
split_go <- function(x) {
  # Taking apart the GO descriptions and GO ID terms
  x <- x %>%
    select(-starts_with("Old")) %>%
    separate(Term, c("GO_term", "GO_ID"), sep = "GO")
  
  # Changing all the floating :#### GO terms to be GO:####
  x <- mutate_if(x, is.character, str_replace_all, pattern = ":", replacement = "GO:")
  
  # Removing the last parentheses in the GO ID 
  x$GO_ID <- x$GO_ID %>% 
    str_replace("\\)", "")
  
  # Removing the last parentheses in the GO term to have clean looking cells 
  x$GO_term <- x$GO_term %>% 
    str_replace("\\($", "")
  return(x)
}

# Read the annotation. Group by gene IDs, and when a gene is duplicated, choose the annotation with the higher bit score. If bit scores are the same, choose the lower E value. If both are the same, leave them be and both annotations can be attached to the differential expression results. Duplicate genes are removed by enrichR anyway. Filter out missing gene name values, rename gene IDs to gene to match the column from the results. Convert gene names to uppercase
annotation <- read_tsv("../gill_annotation_bit50_E1e-6.txt") %>% 
  group_by(`#gene_id`) %>% 
  filter(BitScore == max(BitScore)) %>% 
  filter(Evalue == min(Evalue)) %>% 
  filter(!is.na(gene_name)) %>% 
  rename(gene_id = `#gene_id`) %>% 
  rename(gene = TrinityID) %>% 
  mutate(gene_name = str_to_upper(gene_name))
nrow(distinct(annotation, gene_name))

# Read in the differential expression results, attach the annotation table. 
# pCO2
pco2_unique_genes <- read_tsv("unique_pCO2_genes.txt") %>% 
  left_join(., annotation)
# Count how many annotated genes there are
nrow(distinct(pco2_unique_genes, gene_name))

# Temperature
temp_unique_genes <- read_tsv("unique_temp_genes.txt") %>% 
  left_join(., annotation)
# Count how many annotated genes there are
nrow(distinct(temp_unique_genes, gene_name))

# pCO2 and Temperature
pco2_temp_unique_genes <- read_tsv("unique_pCO2_temp_genes.txt") %>% 
  left_join(., annotation)
# Count how many annotated genes there are
nrow(distinct(pco2_temp_unique_genes, gene_name))


# Run enrichR, searching the biological process, molecular function, and cellular component databases
pCO2_enriched <- enrichr(pco2_unique_genes$gene_name, databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
pCO2_enrich_bio <- as_tibble(pCO2_enriched[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()

# Molecular function
pCO2_enrich_mol <- as_tibble(pCO2_enriched[["GO_Molecular_Function_2021"]])
pCO2_enrich_mol <- dplyr::filter(pCO2_enrich_mol, Adjusted.P.value < 0.05)
print(pCO2_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pCO2_enrich_mol <- split_go(pCO2_enrich_mol)

# Cellular component
pCO2_enrich_cell <- as_tibble(pCO2_enriched[["GO_Cellular_Component_2021"]])
pCO2_enrich_cell <- dplyr::filter(pCO2_enrich_cell, Adjusted.P.value < 0.05)
print(pCO2_enrich_cell$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pCO2_enrich_cell <- split_go(pCO2_enrich_cell)

# Follow the same process of running enrichR and pulling GO terms for the temperature-unique genes
temp_enriched <- enrichr(temp_unique_genes$gene_name, databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
temp_enrich_bio <- as_tibble(temp_enriched[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()

# Molecular function
temp_enrich_mol <- as_tibble(temp_enriched[["GO_Molecular_Function_2021"]])
temp_enrich_mol <- dplyr::filter(temp_enrich_mol, Adjusted.P.value < 0.05)
print(temp_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
temp_enrich_mol <- split_go(temp_enrich_mol)

# Cellular component
temp_enrich_cell <- as_tibble(temp_enriched[["GO_Cellular_Component_2021"]])
temp_enrich_cell <- dplyr::filter(temp_enrich_cell, Adjusted.P.value < 0.05)
print(temp_enrich_cell$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
temp_enrich_cell <- split_go(temp_enrich_cell)

# The same process, for the pCO2 and temperature combined treatment
pco2_temp_enriched <- enrichr(pco2_temp_unique_genes$gene_name, databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
pco2_temp_enrich_bio <- as_tibble(pco2_temp_enriched[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()

# Molecular function
pco2_temp_enrich_mol <- as_tibble(pco2_temp_enriched[["GO_Molecular_Function_2021"]])
pco2_temp_enrich_mol <- dplyr::filter(pco2_temp_enrich_mol, Adjusted.P.value < 0.05)
print(pco2_temp_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pco2_temp_enrich_mol <- split_go(pco2_temp_enrich_mol)

# Cellular component
pco2_temp_enrich_cell <- as_tibble(pco2_temp_enriched[["GO_Cellular_Component_2021"]])
pco2_temp_enrich_cell <- dplyr::filter(pco2_temp_enrich_cell, Adjusted.P.value < 0.05)
print(pco2_temp_enrich_cell$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pco2_temp_enrich_cell <- split_go(pco2_temp_enrich_cell)

# Write out results within each contrast, combined between the 3 GO databases for Revigo. This step is to remove redundant GO terms. Only including the GO ID and adjusted p value, with no column names
# Turn the scientific notation to decimal for input into Revigo
options(scipen=999)
write_delim(as.data.frame(rbind(cbind(pCO2_enrich_bio$GO_ID, pCO2_enrich_bio$Adjusted.P.value), cbind(pCO2_enrich_mol$GO_ID, pCO2_enrich_mol$Adjusted.P.value), cbind(pCO2_enrich_cell$GO_ID, pCO2_enrich_cell$Adjusted.P.value))), "pCO2_GO.txt", delim = "\t", col_names = FALSE)
write_delim(as.data.frame(rbind(cbind(temp_enrich_bio$GO_ID, temp_enrich_bio$Adjusted.P.value), cbind(temp_enrich_mol$GO_ID, temp_enrich_mol$Adjusted.P.value), cbind(temp_enrich_cell$GO_ID, temp_enrich_cell$Adjusted.P.value))), "temp_GO.txt", delim = "\t", col_names = FALSE)
write_delim(as.data.frame(rbind(cbind(pco2_temp_enrich_bio$GO_ID, pco2_temp_enrich_bio$Adjusted.P.value), cbind(pco2_temp_enrich_mol$GO_ID, pco2_temp_enrich_mol$Adjusted.P.value), cbind(pco2_temp_enrich_cell$GO_ID, pco2_temp_enrich_cell$Adjusted.P.value))), "pco2_temp_GO.txt", delim = "\t", col_names = FALSE)


# Write out the full biological process tables for the three treatments for use as supplementary tables in the manuscript
write_tsv(pCO2_enrich_bio, file = "pco2_bio_unique_terms.txt", col_names = TRUE)
write_tsv(temp_enrich_bio, file = "temp_bio_unique_terms.txt", col_names = TRUE)
write_tsv(pco2_temp_enrich_bio, file = "pco2_temp_bio_unique_terms.txt", col_names = TRUE)

# Read in the Revigo results for non-redundant GO terms, only keep class representatives by filtering for null
pco2_revigo_bio <- read_tsv("pCO2_Revigo_BP_Table.tsv") %>% 
  filter(Representative == "null") %>% 
  mutate(database = "Biological Process 2021")
pco2_revigo_mol <- read_tsv("pCO2_Revigo_MF_Table.tsv") %>% 
  filter(Representative == "null") %>% 
  mutate(database = "Molecular Function 2021")
pco2_revigo_cell <- read_tsv("pCO2_Revigo_CC_Table.tsv") %>% 
  filter(Representative == "null") %>% 
  mutate(database = "Cellular Component 2021")

# Make a combined enrichR table
pco2_enrich_combined <- as_tibble(rbind(pCO2_enrich_bio, pCO2_enrich_mol, pCO2_enrich_cell))

# Bind the Revigo results with the enrichR table to retain enrichR p-values and combined scores
# Sort by enrichment database and combined score for plotting
pco2_revigo <- as_tibble(rbind(pco2_revigo_bio, pco2_revigo_mol, pco2_revigo_cell)) %>% 
  left_join(., pco2_enrich_combined, by = c("TermID" = "GO_ID")) %>% 
  arrange(factor(database, levels = c("Cellular Component 2021", "Molecular Function 2021", "Biological Process 2021")), Combined.Score)
# Turn the individual GO term IDs into factors so ggplot respects the order of terms
pco2_revigo$TermID <- factor(pco2_revigo$TermID, levels = pco2_revigo$TermID)

# Repeat the process for the temperature results
temp_revigo_bio <- read_tsv("temp_Revigo_BP_Table.tsv") %>% 
  filter(Representative == "null") %>% 
  mutate(database = "Biological Process 2021")
temp_revigo_mol <- read_tsv("temp_Revigo_MF_Table.tsv") %>% 
  filter(Representative == "null") %>% 
  mutate(database = "Molecular Function 2021")
temp_revigo_cell <- read_tsv("temp_Revigo_CC_Table.tsv") %>% 
  filter(Representative == "null") %>% 
  mutate(database = "Cellular Component 2021")

# Make a combined enrichR table
temp_enrich_combined <- as_tibble(rbind(temp_enrich_bio, temp_enrich_mol, temp_enrich_cell))

# Bind the Revigo results with the enrichR table to retain enrichR p-values and combined scores
# Sort by enrichment database and combined score for plotting
temp_revigo <- as_tibble(rbind(temp_revigo_bio, temp_revigo_mol, temp_revigo_cell)) %>% 
  left_join(., temp_enrich_combined, by = c("TermID" = "GO_ID")) %>% 
  arrange(factor(database, levels = c("Cellular Component 2021", "Molecular Function 2021", "Biological Process 2021")), Combined.Score)
# Turn the individual GO term IDs into factors so ggplot respects the order of terms
temp_revigo$TermID <- factor(temp_revigo$TermID, levels = temp_revigo$TermID)

# Repeat the process for the pco2 and temperature results
pco2_temp_revigo_bio <- read_tsv("pco2_temp_Revigo_BP_Table.tsv") %>% 
  filter(Representative == "null") %>% 
  mutate(database = "Biological Process 2021")
pco2_temp_revigo_mol <- read_tsv("pco2_temp_Revigo_MF_Table.tsv") %>% 
  filter(Representative == "null") %>% 
  mutate(database = "Molecular Function 2021")
pco2_temp_revigo_cell <- read_tsv("pco2_temp_Revigo_CC_Table.tsv") %>% 
  filter(Representative == "null") %>% 
  mutate(database = "Cellular Component 2021")

# Make a combined enrichR table
pco2_temp_enrich_combined <- as_tibble(rbind(pco2_temp_enrich_bio, pco2_temp_enrich_mol, pco2_temp_enrich_cell))

# Bind the Revigo results with the enrichR table to retain enrichR p-values and combined scores
# Sort by enrichment database and combined score for plotting
pco2_temp_revigo <- as_tibble(rbind(pco2_temp_revigo_bio, pco2_temp_revigo_mol, pco2_temp_revigo_cell)) %>% 
  left_join(., pco2_temp_enrich_combined, by = c("TermID" = "GO_ID")) %>% 
  arrange(factor(database, levels = c("Cellular Component 2021", "Molecular Function 2021", "Biological Process 2021")), Combined.Score)
# Turn the individual GO term IDs into factors so ggplot respects the order of terms
pco2_temp_revigo$TermID <- factor(pco2_temp_revigo$TermID, levels = pco2_temp_revigo$TermID)

# Plot the GO terms in a bar plot
pco2_revigo_plot <- ggplot(data = pco2_revigo, aes(x = TermID, y = Combined.Score, colour = database, fill = database, group = database)) +
  geom_col() +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme_classic() +
  xlab(element_blank()) +
  ylab("Combined Score") +
  scale_colour_manual(name = "Enrichment Database", labels = c("Biological Process 2021" = "Biological Process", "Molecular Function 2021" = "Molecular Function", "Cellular Component 2021" = "Cellular Component"), values = c("#4575b4", "#fdae61", "#d73027"), aesthetics = c("colour", "fill"), breaks = c("Biological Process 2021", "Molecular Function 2021", "Cellular Component 2021")) +
  scale_x_discrete(labels = pco2_revigo$Name) +
  theme(text=element_text(size=16), axis.ticks = element_blank()) +
  coord_flip()  

temp_revigo_plot <- ggplot(data = temp_revigo, aes(x = TermID, y = Combined.Score, colour = database, fill = database, group = database)) +
  geom_col() +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme_classic() +
  xlab(element_blank()) +
  ylab("Combined Score") +
  scale_colour_manual(name = "Enrichment Database", labels = c("Biological Process 2021" = "Biological Process", "Molecular Function 2021" = "Molecular Function", "Cellular Component 2021" = "Cellular Component"), values = c("#4575b4", "#fdae61", "#d73027"), aesthetics = c("colour", "fill"), breaks = c("Biological Process 2021", "Molecular Function 2021", "Cellular Component 2021")) +
  scale_x_discrete(labels = temp_revigo$Name) +
  theme(text=element_text(size=16), axis.ticks = element_blank()) +
  coord_flip()  


pco2_temp_revigo_plot <- ggplot(data = pco2_temp_revigo, aes(x = TermID, y = Combined.Score, colour = database, fill = database, group = database)) +
  geom_col() +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme_classic() +
  xlab(element_blank()) +
  ylab("Combined Score") +
  scale_colour_manual(name = "Enrichment Database", labels = c("Biological Process 2021" = "Biological Process", "Molecular Function 2021" = "Molecular Function", "Cellular Component 2021" = "Cellular Component"), values = c("#4575b4", "#fdae61", "#d73027"), aesthetics = c("colour", "fill"), breaks = c("Biological Process 2021", "Molecular Function 2021", "Cellular Component 2021")) +
  scale_x_discrete(labels = pco2_temp_revigo$Name) +
  theme(text=element_text(size=16), axis.ticks = element_blank()) +
  coord_flip()  

ggsave(filename = "pco2_revigo_plot.pdf", plot = pco2_revigo_plot, width = 15, height = 20)
ggsave(filename = "temp_revigo_plot.pdf", plot = temp_revigo_plot, width = 15, height = 23)
ggsave(filename = "pco2_temp_revigo_plot.pdf", plot = pco2_temp_revigo_plot, width = 15, height = 20)

# Write out the three tables 
write_tsv(x = pco2_revigo, file = "pco2_revigo.txt")
write_tsv(x = temp_revigo, file = "temp_revigo.txt")
write_tsv(x = pco2_temp_revigo, file = "pco2_temp_revigo.txt")
