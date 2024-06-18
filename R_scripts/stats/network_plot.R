library(tidyverse) # Manage data
library(fishualize) # Nice plot colours
library(moments) # Calculate skewness
library(enrichR) # Check for GO terms in network results
library(ComplexUpset)
library(patchwork)
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

# Read in the three cluster results, add a column of dataset for combining the data later and an index column for labeling genes as having the highest to lowest betweenness
pCO2_network <- read_csv("pCO2_node_table.csv") %>% 
  mutate(dataset = "pCO2") %>% 
  arrange(desc(Betweenness)) %>% 
  mutate(index = seq(1:n()))

temp_network <- read_csv("temp_node_table.csv") %>% 
  mutate(dataset = "temp") %>% 
  arrange(desc(Betweenness)) %>% 
  mutate(index = seq(1:n()))


combined_network <- read_csv("pCO2_temp_node_table.csv") %>% 
  mutate(dataset = "combo") %>% 
  arrange(desc(Betweenness)) %>% 
  mutate(index = seq(1:n()))


# Combine the datasets, but only with genes with betweenness greater than 0
combined_data <- rbind(temp_network, pCO2_network, combined_network) %>% 
  filter(Betweenness > 0)

# Take a quick look at the data. Most importantly, identify skewness in each dataset
combined_data %>% 
  group_by(dataset) %>% 
  summarise(mean = mean(Betweenness), sd = sd(Betweenness), median = median(Betweenness), n = n(), avg = mean(Betweenness)/n(), cv = mean(Betweenness)/sd(Betweenness), skew = moments::skewness(Betweenness))

# Try out confidence intervals via bootstrapping. They did not work
skewness_fn <- function(data, i){
  d2 <- data[i,]
  return(moments::skewness(d2$Betweenness))
}

skewness_fn(combined_data)


boot::boot.ci(boot::boot(data = filter(combined_data, Betweenness > 0 & dataset == "combo"), statistic = skewness_fn, R = 10000), type = "perc", conf = 0.95)
boot::boot.ci(boot::boot(data = filter(combined_data, Betweenness > 0 & dataset == "temp"), statistic = skewness_fn, R = 10000), type = "perc", conf = 0.95)
boot::boot.ci(boot::boot(data = filter(combined_data, Betweenness > 0 & dataset == "pCO2"), statistic = skewness_fn, R = 10000), type = "perc", conf = 0.95)


# Plot the data with violin and histogram plots
ggplot(data = combined_data, aes(x = dataset, y = Betweenness)) +
  geom_violin() +
  theme_bw()

ggplot(data = combined_data, aes(x = Betweenness, fill = dataset)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(~ dataset)



ggplot(data = pCO2_network, aes(x = Betweenness)) +
  geom_histogram(bins = 50) +
  xlim(-10, 200000) +
  #ylim(0, 300) +
  theme_bw()


ggplot(data = temp_network, aes(x = Betweenness)) +
  geom_histogram(bins = 50) +
  xlim(-10, 200000) +
  #ylim(0, 300) +
  theme_bw()


ggplot(data = combined_network, aes(x = Betweenness)) +
  geom_histogram(bins = 50) +
  xlim(-10, 200000) +
  #ylim(0, 300) +
  theme_bw()



ggplot(data = temp_network %>% arrange(desc(Betweenness)) %>% filter(Betweenness > 0), aes(x = seq(1:nrow(temp_network %>% filter(Betweenness > 0))), y = Betweenness)) +
  geom_point() +
  theme_bw()
  
# Try separate plots of betweenness and gene index
ggplot(data = pCO2_network %>% arrange(desc(Betweenness)) %>% filter(Betweenness > 0), aes(x = seq(1:nrow(pCO2_network %>% filter(Betweenness > 0))), y = Betweenness)) +
  geom_point() +
  theme_bw()

ggplot(data = combined_network %>% arrange(desc(Betweenness)) %>% filter(Betweenness > 0), aes(x = seq(1:nrow(combined_network %>% filter(Betweenness > 0))), y = Betweenness)) +
  geom_point() +
  theme_bw()

# Make a combined skewness plot, inputting the skewness values identified above
skewness_plot <- ggplot(data = combined_data %>% arrange(dataset, desc(Betweenness)) %>% filter(Betweenness > 0, index < 250), aes(x = index, y = Betweenness, color = dataset, fill = dataset)) +
  geom_point(alpha = 0.9, size= 1.1) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = 1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp*", 3.40"), bquote(italic(p)*CO[2]*", 4.22"), "Temperature, 5.02")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = 1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp*", 3.40"), bquote(italic(p)*CO[2]*", 4.22"), "Temperature, 5.02")) +
  xlab("Gene Index") +
  guides(fill = guide_legend(title="Treatment, Skewness"), colour = guide_legend(title="Treatment, Skewness")) +
  theme_bw(base_size = 18) +
  theme(legend.position = c(0.80, 0.75), legend.background = element_blank(), legend.key = element_blank())
skewness_plot
ggsave(filename = "skew_plot.pdf", plot = skewness_plot, dpi = 2000)

#####################################################################################################################
# Try using enrichR to identify GO terms in the different datasets
# Run enrichR, searching the biological process database
# Run temperature first. First check genes with betweenness of 0 (the 'spokes')
temp_enriched_spokes <- enrichr(temp_network %>% filter(Betweenness == 0) %>% pull(Label), databases = c("GO_Biological_Process_2021"))
# Filter for significant GO terms, sort by combined score, use my custom function to split GO term and GO ID
temp_enriched_spokes_bio <- as_tibble(temp_enriched_spokes[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()

# Then temperature with betweenness > 0 (the 'hubs')
temp_enriched_hubs <- enrichr(temp_network %>% filter(Betweenness > 0) %>% pull(Label), databases = c("GO_Biological_Process_2021"))
# Filter for significant GO terms, sort by combined score, use my custom function to split GO term and GO ID
temp_enriched_hubs_bio <- as_tibble(temp_enriched_hubs[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()

# Now the pCO2 genes. First the spokes, then the hubs
pCO2_enriched_spokes <- enrichr(pCO2_network %>% filter(Betweenness == 0) %>% pull(Label), databases = c("GO_Biological_Process_2021"))
# Filter for significant GO terms, sort by combined score, use my custom function to split GO term and GO ID
pCO2_enrich_spokes_bio <- as_tibble(pCO2_enriched_spokes[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()

pCO2_enriched_hubs <- enrichr(pCO2_network %>% filter(Betweenness > 0) %>% pull(Label), databases = c("GO_Biological_Process_2021"))
# Filter for significant GO terms, sort by combined score, use my custom function to split GO term and GO ID
pCO2_enrich_hubs_bio <- as_tibble(pCO2_enriched_hubs[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()

# Now the pCO2+temp genes. Spokes then hubs
combined_enriched_spokes <- enrichr(combined_network %>% filter(Betweenness == 0) %>% pull(Label), databases = c("GO_Biological_Process_2021"))
# Filter for significant GO terms, sort by combined score, use my custom function to split GO term and GO ID
combined_enrich_spokes_bio <- as_tibble(combined_enriched_spokes[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()

combined_enriched_hubs <- enrichr(combined_network %>% filter(Betweenness > 0) %>% pull(Label), databases = c("GO_Biological_Process_2021"))
# Filter for significant GO terms, sort by combined score, use my custom function to split GO term and GO ID
combined_enrich_hubs_bio <- as_tibble(combined_enriched_hubs[["GO_Biological_Process_2021"]]) %>% 
  dplyr::filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score)) %>% 
  split_go()


# Write out these network results as supplementary tables
write_tsv(x = temp_enriched_spokes_bio, file = "Network_Temp_spokes.txt")
write_tsv(x = temp_enriched_hubs_bio, file = "Network_Temp_hubs.txt")

write_tsv(x = pCO2_enrich_spokes_bio, file = "Network_pCO2_spokes.txt")
write_tsv(x = pCO2_enrich_hubs_bio, file = "Network_pCO2_hubs.txt")

write_tsv(x = combined_enrich_spokes_bio, file = "Network_Combined_spokes.txt")
write_tsv(x = combined_enrich_hubs_bio, file = "Network_Combined_hubs.txt")


# Make an upset plot
# Create an upset plot of the combined and unique term descriptions from the KEGG analysis
combined_list <- list(Temp = c(combined_data %>% filter(dataset == "temp") %>% column_to_rownames(Label) %>%pull(Label)), pCO2 = c(combined_data %>% filter(dataset == "pCO2") %>% pull(Label)), combined = c(combined_data %>% filter(dataset == "combo") %>% pull(Label)))

# Make a complex upset plot rather than the original version
upset_plot <- upset(data = UpSetR::fromList(combined_list), intersect = c("Temp", "pCO2", "combined"), labeller=ggplot2::as_labeller(c('Temp' = 'Temperature','pCO2' = "pCO2", 'combined' =  'Combined')), name = bquote("Treatment"), set_sizes=(upset_set_size()+ylab("Number of Transcripts")))
upset_plot


combined_data %>% distinct(Label) %>% nrow(.)


# Create a dataframe with presence absense information not using UpSetR's function, so we can add extra information to it. First, calculate betweenness centrality scores for each gene based on mean values if there are multiple, or just one if not
melted_data <- combined_data %>% 
  add_column(Present=1) %>% 
  pivot_wider(names_from=dataset, values_from=Present) %>% 
  replace(is.na(.), 0) %>% 
  group_by(Label) %>% 
  summarize(mean_betweenness = mean(Betweenness))

# Create the presence absence dataframe of genes
melted_data2 <- combined_data %>% 
  add_column(Present=1) %>% 
  pivot_wider(names_from=dataset, values_from=Present) %>% 
  replace(is.na(.), 0) %>% 
  group_by(Label) %>% 
  summarize(pCO2 = sum(pCO2), temp = sum(temp), combo = sum(combo))

# Join the mean betweenness centrality scores to the presence abscense data
presence_absense_data <- left_join(melted_data, melted_data2, by = "Label")

# Check that the basic upset plot looks the same- it does.
upset(data = presence_absense_data, intersect = c("temp", "pCO2", "combo"), labeller=ggplot2::as_labeller(c('temp' = 'Temperature','pCO2' = "pCO2", 'combo' =  'Combined')), name = bquote("Treatment"), set_sizes=(upset_set_size()+ylab("Number of Transcripts")), annotations = list('mean_betweenness'=upset_annotate('mean_betweenness', geom_violin(na.rm=TRUE))+ylab(bquote(atop("Mean Betweenness", "Centrality Score")))))

t.test(combined_data %>% filter(dataset == "temp") %>% pull(Betweenness), combined_data %>% filter(dataset == "combo") %>% pull(Betweenness))

ggplot(combined_data, aes(x = dataset, y = Betweenness)) +
  #geom_point() +
  geom_boxplot() +
  theme_bw()


upset_plot_networktranscripts <- upset(data = presence_absense_data, intersect = c("temp", "pCO2", "combo"), labeller=ggplot2::as_labeller(c('temp' = 'Temperature','pCO2' = "pCO2", 'combo' =  'Combined')), name = bquote("Treatment"), set_sizes=(upset_set_size()+ylab("Number of Transcripts")))
ggsave(filename = "upsetplot_network.pdf", plot = upset_plot_networktranscripts, dpi = 3000)
ggsave(filename = "upsetplot_network.svg", plot = upset_plot_networktranscripts, dpi = 3000)


skew_upset_plot <- skewness_plot / upset_plot_networktranscripts
ggsave(filename = "upsetplot_skew_network.svg", plot = skew_upset_plot, dpi = 3000)
