library(tidyverse)
library(brms) # Run the Bayesian models
library(tidybayes) # Plot Bayesian models
library(fishualize)
library(patchwork)
# Faster Hamiltonian Monte Carlo sampling
options(brms.backend = "cmdstanr")

# Read in data following the data input tutorial https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf
# Read the log2-transformed counts per million data from EdgeR
read_data <- read_tsv("../edgeR_DGE/cpm_nolog.txt") %>% 
  column_to_rownames(var = "Gene_ID")
dim(read_data)

# Create a vector of NKA transripts
nka_transcripts <- unique(c("TR55914|c0_g1_i1",
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
                     "TR70983|c0_g1_i3"))
nka_counts <- as_tibble(read_data[nka_transcripts,], rownames = NA) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_ID")

# Read in the ATPase data. Take the mean of the 2 ATPase activity measurements for each individual
atpase_mean <- readxl::read_excel("../Physiology/NaK-ATPase Data Update Sept 1 2022.xlsx", sheet = "NaK-ATPase_reformat") %>% 
  separate(Treatment, into = c("Treatment", "Index", "Time")) %>% 
  mutate(Time = case_when(Time == "Recovery" ~ 720,
                          Time == "Baseline" ~ 0, 
                          Time == "24" ~ 24,
                          Time == "7" ~ 168,
                          Time == "6" ~ 6)) %>% 
  rename(ATPase = `nmol NADH/min*mgpr`, Run = `Run #`) %>% 
  mutate(treatment_combined = paste0(Treatment, "_", Time), Ind_ID = paste0(Treatment, "_", Time, "_", Index)) %>% 
  filter(!is.na(RNA_ID)) %>% 
  group_by(RNA_ID) %>% 
  summarize(Mean_ATPase = mean(ATPase))
  
# Read in the ATPase data again to get the experiment metadata 
atpase <- readxl::read_excel("../Physiology/NaK-ATPase Data Update Sept 1 2022.xlsx", sheet = "NaK-ATPase_reformat") %>% 
  separate(Treatment, into = c("Treatment", "Index", "Time")) %>% 
  mutate(Time = case_when(Time == "Recovery" ~ 720,
                          Time == "Baseline" ~ 0, 
                          Time == "24" ~ 24,
                          Time == "7" ~ 168,
                          Time == "6" ~ 6)) %>% 
  rename(ATPase = `nmol NADH/min*mgpr`, Run = `Run #`) %>% 
  mutate(treatment_combined = paste0(Treatment, "_", Time), Ind_ID = paste0(Treatment, "_", Time, "_", Index)) %>% 
  filter(!is.na(RNA_ID)) %>% 
  select(-'Δ AOuabain', -'ΔA ATPase', -'Net Δ Abs', -'Protein (mg/mL)', -'ATPase') %>% 
  distinct(RNA_ID, .keep_all = TRUE)

# Join the metadata to the mean values
atpase_mean <- left_join(atpase_mean, atpase, by = "RNA_ID")

# Reorder the atpase data to match the read data. Discard individuals in the ATPase data missing from the read data
nka_atpase <- left_join(atpase_mean, nka_counts, by = "RNA_ID") %>% 
  filter(complete.cases(.))

# Clean up the environment now that we have a combined dataframe
rm(atpase, atpase_mean, nka_counts, read_data, nka_transcripts)

# Replace the | character in column names with . to keep brms happy
names(nka_atpase) <- gsub(x = names(nka_atpase), pattern = "\\|", replacement = ".")  

# # Test a linear model 
# tr1_model <- brm(data = nka_atpase, family = gaussian,
#                  Mean_ATPase ~ TR74551.c2_g1_i4, 
#                           prior = c(prior(normal(0, 80), class = b),
#                                     prior(normal(0, 80), class = Intercept)),
#                           iter = 5000, warmup = 1000, chains = 4, cores = 4,  
#                           control = list(adapt_delta = 0.98, max_treedepth = 12))
# 
# summary(tr1_model)
# bayes_R2(tr1_model)
# plot(conditional_effects(tr1_model), points = T, theme = theme_bw())

# Take only the transcript IDs from the NKA data table
transcript_IDs <- colnames(nka_atpase)[12:44]

# Create a function that runs a Bayesian model, summarizes it, and identifies R2, and outputs a list
model_function <- function(gene, dataset){
  model <- brm(data = dataset, family = student,
                   bf(paste("Mean_ATPase ~ ", gene)), 
                   prior = c(prior(normal(0, 100), class = b),
                             prior(normal(0, 100), class = Intercept)),
                   iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                   control = list(adapt_delta = 0.98, max_treedepth = 12), 
                   save_pars = save_pars(all = TRUE))
  model_summary <- summary(model)
  R2 <- bayes_R2(model)
  model <- add_criterion(model, "loo", moment_match = TRUE)
  model_loo <- loo(model, moment_match = TRUE)
  model_list <- list(model = model, summary = model_summary, R2 = R2, model_loo = model_loo)
  return(model_list)
}


# Initialize an empty list
NKA_transcript_models <- list()

# Loop through NKA transcripts and create linear models for each. Add results to the empty overall list
for (i in 1:length(transcript_IDs)){
  
  print(transcript_IDs[i])# Print the gene being worked on
  
  NKA_transcript_models[[transcript_IDs[i]]] <- model_function(transcript_IDs[i], nka_atpase) # Use the model function with the present transcript
  
}

rm(i)

# Create an empty object, then loop through the list of models and pull estimates for correlation coefficients and coefficients of determination
model_summaries_combined <- c()

for (i in 1:length(transcript_IDs)){
  
  print(transcript_IDs[i])# Print the gene being worked on
  
  transcript_data <- NKA_transcript_models[[i]]
  
  model_intercept <- transcript_data$summary$fixed[1, ]
  colnames(model_intercept) <- c("Intercept.Estimate", "Intercept.Error", "Intercept.l-95% CI", "Intercept.u-95% CI", "Int.Rhat", "Int.Bulk_ESS", "Int.Tail_ESS")
  
  model_summary <- transcript_data$summary$fixed[2, ]
  colnames(model_summary) <- c("Transcript.Estimate", "Transcript.Error", "Transcript.l-95% CI", "Transcript.u-95% CI", "Transcript.Rhat", "Transcript.Bulk_ESS", "Transcript.Tail_ESS")
  
  model_R2 <- transcript_data$R2
  colnames(model_R2) <- c("R2", "R2.Error", "R2.Q2.5", "R2.Q97.5")
  
  
  
  model_elpd <- transcript_data$model_loo$estimates[1,1]
  model_elpd.se <- transcript_data$model_loo$estimates[1,2]
  model_loo_summ <- cbind(model_elpd, model_elpd.se)
  colnames(model_loo_summ) <- c("elpd_loo", "loo.se")
  
  model_out <- cbind(model_intercept, model_summary, model_R2, model_loo_summ)
  
  rownames(model_out) <- rownames(model_summary)
  
  model_summaries_combined <- rbind(model_summaries_combined, model_out)
}
rm(i, transcript_data, model_intercept, model_summary, model_R2, model_out, model_elpd, model_elpd.se, model_loo_summ)

# Sort by descending transcript estimates and R2 values
model_summaries_combined <- model_summaries_combined %>% 
  arrange(desc(elpd_loo), desc(R2), desc(Transcript.Estimate)) %>% 
  rownames_to_column(var = "Transcript_ID")

# Replace the , character in column names with | to be consistent with annotations
model_summaries_combined$Transcript_ID <- gsub(x = model_summaries_combined$Transcript_ID, pattern = "\\.", replacement = "|") 


# Read the annotation. Group by gene IDs, and when a gene is duplicated, choose the annotation with the higher bit score. If bit scores are the same, choose the lower E value. If both are the same, leave them be and both annotations can be attached to the differential expression results. Duplicate genes are removed by enrichR anyway. Filter out missing gene name values, rename gene IDs to gene to match the column from the results. Convert gene names to uppercase
annotation <- read_tsv("../gill_annotation_bit50_E1e-6.txt") %>% 
  filter(!is.na(sp.BX_transcript_name)) %>%
  group_by(TrinityID) %>% 
  arrange(desc(BitScore), Evalue) %>% 
  slice_head(n = 1) %>% 
  rename(Transcript_ID = TrinityID, gene_id = `#gene_id`) %>% 
  mutate(gene_name = str_to_upper(gene_name)) %>%
  select(Transcript_ID, gene_id, gene_name, sp.BX_transcript_name)



model_annotations_combined <- left_join(model_summaries_combined, annotation, by = "Transcript_ID")


plot(conditional_effects(NKA_transcript_models$TR74551.c2_g1_i9$model), points = T, theme = theme_bw(base_size = 14))

plot1 <- plot(conditional_effects(NKA_transcript_models$TR70983.c0_g1_i1$model), points = F, theme = theme_bw(base_size = 14))

best_model_plot <- plot1$TR70983.c0_g1 +
  xlab(bquote(atop("Subunit "*italic(beta)*"3 (TR70983|c0_g1_i1)", "Counts per Million"))) +
  ylab(bquote(atop(Na^"+"*"/"*K^"+"~ATPase~Activity, (nmol~NADH*"\U00B7"*min^-1*"\U00B7"*mg~pr^-1)))) +
  geom_point(data = nka_atpase, aes(x = `TR70983.c0_g1_i1`, y = Mean_ATPase, colour = Treatment, shape = factor(Time)), inherit.aes = FALSE) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control"), breaks = c("Combo", "pCO2", "Temp", "Control")) +
  guides(colour = "none", shape = guide_legend("Timepoint")) +
  theme(legend.background = element_blank(), legend.position = c(.86, 0.76))
best_model_plot

plot2 <- plot(conditional_effects(NKA_transcript_models$TR45365.c0_g1_i3$model), points = T, theme = theme_bw(base_size = 18))
plot2$TR45365.c0_g1_i3 +
  xlab("TR45365|c0_g1_i3 Counts per Million") +
  ylab(bquote(Na^"+"*"/"*K^"+"~ATPase~Activity~(nmol~NADH*"\U00B7"*min^-1*"\U00B7"*mg~pr^-1)))
  

plot(conditional_effects(NKA_transcript_models$TR85165.c0_g1_i1$model), points = T, theme = theme_bw())

# Add a new column with shortened subunit names for nicer plotting
model_annotations_combined <- model_annotations_combined %>% 
  separate(sp.BX_transcript_name, sep = " ", into = c("Sodium", "ATPase", "foo", "subunit"), remove = FALSE) %>% 
  select(-Sodium, -ATPase, -foo)  %>% 
  separate(subunit, into = c("subunit_overall", "number"), sep = "-", remove = FALSE) %>% 
  select(-number)


ggplot(model_annotations_combined, aes(y = Transcript.Estimate, x = subunit, ymin = `Transcript.l-95% CI`, ymax = `Transcript.u-95% CI`)) +
  geom_pointrange(position=position_dodge2(width=0.70)) +
  ylim(-50, 50) +
  theme_bw()


ggplot(model_annotations_combined, aes(y = Transcript.Estimate, x = elpd_loo, color = subunit)) +
  geom_point() +
  theme_bw()



subunit_R2_plot <- ggplot(model_annotations_combined, aes(x = subunit, y = R2)) +
  geom_violin() +
  geom_point() +
  #geom_errorbar(aes(ymin = R2.Q2.5, ymax = R2.Q97.5)) +
  ylab(bquote(italic(R)^2)) +
  xlab(bquote(Na^"+"*"/"*K^"+"~ATPase~Subunit)) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
subunit_R2_plot
ggsave(filename = "NKA_subunit_R2.pdf", plot = subunit_R2_plot, dpi = 2000)

# Plot the loo values for each model, separated by subunit
subunit_loo_plot <- ggplot(model_annotations_combined, aes(x = subunit, y = elpd_loo)) +
  geom_violin() +
  geom_point() +
  #geom_errorbar(aes(ymin = elpd_loo - loo.se, ymax = elpd_loo + loo.se)) +
  ylab("ELPD Loo") +
  xlab(bquote(Na^"+"*"/"*K^"+"~ATPase~Subunit)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
subunit_loo_plot
ggsave(filename = "NKA_subunit_loo.pdf", plot = subunit_loo_plot, dpi = 2000)


subunit_loo_plot_simplified <- ggplot(model_annotations_combined, aes(x = subunit_overall, y = elpd_loo)) +
  geom_boxplot(outliers = FALSE) +
  ggbeeswarm::geom_beeswarm() +
  #geom_errorbar(aes(ymin = elpd_loo - loo.se, ymax = elpd_loo + loo.se)) +
  ylab("ELPD LOO") +
  xlab(bquote(Na^"+"*"/"*K^"+"~ATPase~Subunit)) +
  scale_x_discrete(labels = c(bquote(alpha), bquote(beta))) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
subunit_loo_plot_simplified

model_annotations_combined %>% filter(subunit_overall == "alpha") %>% pull(elpd_loo) %>% mean(.)
model_annotations_combined %>% filter(subunit_overall == "beta") %>% pull(elpd_loo) %>% mean(.)

model_annotations_combined %>% filter(subunit_overall == "alpha") %>% pull(elpd_loo) %>% mean(.)
model_annotations_combined %>% filter(subunit_overall == "beta") %>% pull(elpd_loo) %>% mean(.)

model_annotations_combined %>% 
  group_by(subunit_overall) %>% 
  summarize(mean_loo = mean(Transcript.Estimate))


write_tsv(x = model_annotations_combined, file = "NKA_mRNA_activity_models.txt")

#####################################################################################################################################################
# Following a suggestion by Trish Schulte, run the different models with a variable selection approach to find the model with the most predictive power for NKA activity. Models will be run by descending values of loo above


model_1var <- brm(data = nka_atpase, family = gaussian,
             bf(Mean_ATPase ~ TR74551.c2_g1_i3), 
             prior = c(prior(normal(0, 100), class = b),
                       prior(normal(0, 100), class = Intercept)),
             iter = 20000, warmup = 5000, chains = 4, cores = 4,  
             control = list(adapt_delta = 0.98, max_treedepth = 12), 
             save_pars = save_pars(all = TRUE))
model_1var <- add_criterion(model_1var, "loo", moment_match = TRUE)
summary(model_1var)


model_2var <- brm(data = nka_atpase, family = gaussian,
                  bf(Mean_ATPase ~ TR74551.c2_g1_i3 + TR85165.c0_g1_i1), 
                  prior = c(prior(normal(0, 100), class = b),
                            prior(normal(0, 100), class = Intercept)),
                  iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                  control = list(adapt_delta = 0.98, max_treedepth = 12), 
                  save_pars = save_pars(all = TRUE))
model_2var <- add_criterion(model_2var, "loo", moment_match = TRUE)
summary(model_2var)

model_3var <- brm(data = nka_atpase, family = gaussian,
                  bf(Mean_ATPase ~ TR74551.c2_g1_i3 + TR85165.c0_g1_i1 + TR74551.c2_g1_i1), 
                  prior = c(prior(normal(0, 100), class = b),
                            prior(normal(0, 100), class = Intercept)),
                  iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                  control = list(adapt_delta = 0.98, max_treedepth = 12), 
                  save_pars = save_pars(all = TRUE))
model_3var <- add_criterion(model_3var, "loo", moment_match = TRUE)
summary(model_3var)


model_4var <- brm(data = nka_atpase, family = gaussian,
                  bf(Mean_ATPase ~ TR74551.c2_g1_i3 + TR85165.c0_g1_i1 + TR74551.c2_g1_i1 + TR22292.c1_g1_i1), 
                  prior = c(prior(normal(0, 100), class = b),
                            prior(normal(0, 100), class = Intercept)),
                  iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                  control = list(adapt_delta = 0.98, max_treedepth = 12), 
                  save_pars = save_pars(all = TRUE))
model_4var <- add_criterion(model_4var, "loo", moment_match = TRUE)
summary(model_4var)

loo_compare(model_1var, model_2var, model_3var, model_4var)

loo_list <- list(model_1var$criteria$loo, model_2var$criteria$loo, model_3var$criteria$loo, model_4var$criteria$loo)
loo_model_weights(loo_list,method = c("stacking"))

# Load the ATPase plot from the models of ATPase activity over time
load(file = "atpase_plot.RData")

atpase_plot <- ggplot(data = atpase_emmeansdraws, aes(x = .value, y = time, fill = treatment, colour = treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.5) +
  #geom_point(data = raw_data, aes(x = ATPase, y = Time, fill = Treatment, colour = Treatment), alpha = 0.5, shape = 15) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control"), breaks = c("Combo", "pCO2", "Temp", "Control")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control"), breaks = c("Combo", "pCO2", "Temp", "Control")) +
  xlab(bquote(Na^"+"*"/"*K^"+"~ATPase~Activity~(nmol~NADH*"\U00B7"*min^-1*"\U00B7"*mg~pr^-1))) +
  ylab("Timepoint (hours)") +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment")) +
  xlim(4, 42) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.92, 0.30), legend.background = element_blank(), legend.text=element_text(size=12))
atpase_plot
# Make a combined plot of the different ATPase activity and transcript models. The activity over time and ELPD Loo plots are likely the most important, while the particular subunit plot was chosen for its lowest LOO and highest R^2 value
combined_atpase <- atpase_plot / (best_model_plot | subunit_loo_plot_simplified) + plot_annotation(tag_levels = 'A')

ggsave(filename = "combined_ATPase_models.pdf", plot = combined_atpase, dpi = 3000, scale = 2)


###################################################################
# Test the TR74551|c2_g1_i3 model without the final 2 datapoints
# The student distribution was more robust to outliers
subunit_B1_outlier <- brm(data = nka_atpase, family = student,
             bf(Mean_ATPase ~ TR74551.c2_g1_i3), 
             prior = c(prior(normal(0, 100), class = b),
                       prior(normal(0, 100), class = Intercept)),
             iter = 20000, warmup = 5000, chains = 4, cores = 4,  
             control = list(adapt_delta = 0.98, max_treedepth = 12), 
             save_pars = save_pars(all = TRUE))

summary(subunit_B1_outlier)
plot(conditional_effects(subunit_B1_outlier), points = T, theme = theme_bw(base_size = 18)) +
  ylab(bquote(Na^"+"*"/"*K^"+"~ATPase~Activity~(nmol~NADH*"\U00B7"*min^-1*"\U00B7"*mg~pr^-1))) +
  xlab("TR45365|c0_g1_i3 Counts per Million")

