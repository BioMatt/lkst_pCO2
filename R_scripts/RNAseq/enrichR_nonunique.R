library(tidyverse)
library(enrichR)
library(ggdist)

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
  mutate(gene_name_original = gene_name) %>% 
  mutate(gene_name = str_to_upper(gene_name)) 
nrow(distinct(annotation, gene_name))

# Read in the differential expression results, attach the annotation table. 
# pCO2
pco2_genes <- read_tsv("sig_pCO2_DGE_nonunique.txt") %>% 
  left_join(., annotation)
# Count how many annotated genes there are
nrow(distinct(pco2_genes, gene_name))

# Filter for 7 days vs 0 hours
pco2_genes_168v0 <- filter(pco2_genes, contrast == "pCO2_168v0")
# Filter for 6 hours vs 0 hours
pco2_genes_6v0 <- filter(pco2_genes, contrast == "pCO2_6v0")


# Temperature
temp_genes <- read_tsv("sig_temp_DGE_nonunique.txt") %>% 
  left_join(., annotation)
# Count how many annotated genes there are
nrow(distinct(temp_genes, gene_name))

# pCO2 and Temperature
pco2_temp_genes <- read_tsv("sig_pCO2_temp_DGE_nonunique.txt") %>% 
  left_join(., annotation)
# Count how many annotated genes there are
nrow(distinct(pco2_temp_genes, gene_name))

# Run enrichR, searching the biological process, molecular function, and cellular component databases
pCO2_enriched <- enrichr(pco2_genes$gene_name, databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pCO2_enrich_bio <- as_tibble(pCO2_enriched[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()
print(pCO2_enrich_bio$Term)
# plus-end-directed vesicle transport along microtubule (GO:0072383) had a ridiculously high combined score compared to all others

# Run enrichR, searching the biological process, molecular function, and cellular component databases
pCO2_enriched_168v0 <- enrichr(pco2_genes_168v0$gene_name, databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
pCO2_168v0_enrich_bio <- as_tibble(pCO2_enriched_168v0[["GO_Biological_Process_2021"]])
pCO2_168v0_enrich_bio <- dplyr::filter(pCO2_168v0_enrich_bio, Adjusted.P.value < 0.05)
print(pCO2_168v0_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pCO2_168v0_enrich_bio <- split_go(pCO2_168v0_enrich_bio)

# Run enrichR, searching the biological process, molecular function, and cellular component databases
pCO2_enriched_6v0 <- enrichr(pco2_genes_6v0$gene_name, databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
pCO2_6v0_enrich_bio <- as_tibble(pCO2_enriched_6v0[["GO_Biological_Process_2021"]])
pCO2_6v0_enrich_bio <- dplyr::filter(pCO2_6v0_enrich_bio, Adjusted.P.value < 0.05)
print(pCO2_6v0_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pCO2_6v0_enrich_bio <- split_go(pCO2_6v0_enrich_bio)


# Run enrichR, searching the biological process, molecular function, and cellular component databases
temp_enriched <- enrichr(temp_genes$gene_name, databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
temp_enrich_bio <- as_tibble(temp_enriched[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()
print(temp_enrich_bio$Term)

# Run enrichR, searching the biological process, molecular function, and cellular component databases
pCO2_temp_enriched <- enrichr(pco2_temp_genes$gene_name, databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pCO2_temp_enrich_bio <- as_tibble(pCO2_temp_enriched[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()
print(pCO2_temp_enrich_bio$Term)
# Negative regulation of organ growth (GO:0046621) had a ridiculously high combined score, especially compared to the next highest

# Write out one gene list for each of the three treatments for Will to look at the gene network stuff with
write_tsv(as_tibble(pco2_temp_genes$gene_name_original) %>% filter(!is.na(value)) %>% distinct(), file = "pco2_temp_genes.txt", col_names = FALSE)
write_tsv(as_tibble(pco2_genes$gene_name_original) %>% filter(!is.na(value)) %>% distinct(), file = "pco2_genes.txt", col_names = FALSE)
write_tsv(as_tibble(temp_genes$gene_name_original) %>% filter(!is.na(value)) %>% distinct(), file = "temp_genes.txt", col_names = FALSE)

# Write out the different GO term results in tables
write_tsv(pCO2_enrich_bio, file = "pco2_bio_terms.txt", col_names = TRUE)
write_tsv(temp_enrich_bio, file = "temp_bio_terms.txt", col_names = TRUE)
write_tsv(pCO2_temp_enrich_bio, file = "pco2_temp_bio_terms.txt", col_names = TRUE)

# A function to take a table of GO term results from enrichR, and look up the log fold change values from a contrast table from edgeR
GO_foldchange <- function(enrichR_table, contrast_table){
  gene_FCs <- tibble(GO_term = character(), GO_ID = character(), gene_name = character(), logFC = numeric(), abs_logFC = numeric())
  for (GO in 1:nrow(enrichR_table)){
    genes <- enrichR_table[GO, 8]
    genes <- str_split(string = genes, pattern = "\\;")
    genes <- as.vector(genes[[1]])
    for (gene in 1:length(genes)){
      gene_row <- contrast_table[grepl(paste0("^", genes[gene]), contrast_table$gene_name),]
      gene_FCs <- gene_FCs %>% add_row(GO_term = rep(as.character(enrichR_table[GO, 1]), nrow(gene_row)), GO_ID = rep(as.character(enrichR_table[GO, 2]), nrow(gene_row)), gene_name = rep(as.character(genes[gene]), nrow(gene_row)), logFC = gene_row$logFC, abs_logFC = abs(gene_row$logFC))
    }
  }
return(gene_FCs)
}

# Use the custom function to attach fold changes to the pCO2 168 v 0 results
pCO2_168v0_GO_FC <- GO_foldchange(enrichR_table = pCO2_168v0_enrich_bio, contrast_table = pco2_genes_168v0)

pCO2_168v0_GOgeneFCs <- left_join(pCO2_168v0_enrich_bio, pCO2_168v0_GO_FC)
pCO2_168v0_GOgeneFCs$GO_term <- as.factor(pCO2_168v0_GOgeneFCs$GO_term)
pCO2_168v0_GOgeneFCs$GO_ID <- as.factor(pCO2_168v0_GOgeneFCs$GO_ID)

pCO2_168v0_ion_GOfoldchanges <- pCO2_168v0_GOgeneFCs[grepl(" ion ", pCO2_168v0_GOgeneFCs$GO_term),]

ggplot(data = pCO2_168v0_ion_GOfoldchanges, aes(x = abs_logFC, y = GO_term, group = GO_term)) +
  stat_halfeye() +
  labs(x = "Absolute log-Fold Change", y = "Gene Ontology Term") +
  theme_bw()

pCO2_168v0_ion_GO <- pCO2_168v0_enrich_bio[grepl(" ion ", pCO2_168v0_enrich_bio$GO_term),]


ggplot(data = pCO2_168v0_ion_GO, aes(x = Combined.Score, y = GO_term, group = GO_term)) +
  geom_col() +
  labs(x = "Combined Score", y = element_blank()) +
  theme_bw()

# Use the custom function to attach fold changes to the pCO2 6 v 0 results
pCO2_6v0_GO_FC <- GO_foldchange(enrichR_table = pCO2_6v0_enrich_bio, contrast_table = pco2_genes_6v0)

pCO2_6v0_GOgeneFCs <- left_join(pCO2_6v0_enrich_bio, pCO2_6v0_GO_FC)
pCO2_6v0_GOgeneFCs$GO_term <- as.factor(pCO2_6v0_GOgeneFCs$GO_term)
pCO2_6v0_GOgeneFCs$GO_ID <- as.factor(pCO2_6v0_GOgeneFCs$GO_ID)

pCO2_6v0_ion_GOfoldchanges <- pCO2_6v0_GOgeneFCs[grepl(" ion ", pCO2_6v0_GOgeneFCs$GO_term),]

ggplot(data = pCO2_6v0_ion_GOfoldchanges, aes(x = abs_logFC, y = GO_term, group = GO_term)) +
  stat_halfeye() +
  labs(x = "Absolute log-Fold Change", y = "Gene Ontology Term") +
  theme_bw()

pCO2_6v0_ion_GO <- pCO2_6v0_enrich_bio[grepl(" ion ", pCO2_6v0_enrich_bio$GO_term),]


ggplot(data = pCO2_6v0_ion_GO, aes(x = Combined.Score, y = GO_term, group = GO_term)) +
  geom_col() +
  labs(x = "Combined Score", y = element_blank()) +
  theme_bw()

