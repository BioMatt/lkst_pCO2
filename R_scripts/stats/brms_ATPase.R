# Rstan needed to be installed with this line to make the package work
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))b

library(tidyverse) # Handling large datasets
library(brms) # Run the Bayesian models
library(tidybayes) # Plot Bayesian models
library(modelr) # Make model functions work in pipes
library(fishualize) # Get nice fish-based colours
library(DHARMa) # Check residuals and dispersion for Bayesian models
library(performance) # Check model performance
library(modelbased) # Visualize the model's posterior predictions
library(rcartocolor) # Some colours specific to plotting the WAIC following advice from this link: https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/overfitting-regularization-and-information-criteria.html#using-information-criteria
library(bayestestR) # Gives access to the describe_posterior() function
library(emmeans) # Comparisons among groups

# Faster Hamiltonian Monte Carlo sampling
options(brms.backend = "cmdstanr")

# Read in a different version of the raw data, reformatted into long format
raw_data <- readxl::read_excel("../NaK-ATPase Data Update Sept 1 2022.xlsx", sheet = "NaK-ATPase_reformat") %>% 
  separate(Treatment, into = c("Treatment", "Index", "Time")) %>% 
  mutate(time_factor = Time) %>% 
  mutate(Time = case_when(Time == "Recovery" ~ 720,
                          Time == "Baseline" ~ 0, 
                          Time == "24" ~ 24,
                          Time == "7" ~ 168,
                          Time == "6" ~ 6)) %>% 
  rename(ATPase = `nmol NADH/min*mgpr`, Run = `Run #`) %>% 
  mutate(treatment_combined = paste0(Treatment, "_", Time), Ind_ID = paste0(Treatment, "_", Time, "_", Index)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Combo", "pCO2", "Temp", "Control"))) |> 
  mutate(Time = factor(Time, levels = c("0", "6", "24", "168", "720"))) 

raw_data$Run <- as.factor(raw_data$Run)
raw_data$treatment_combined <- as.factor(raw_data$treatment_combined)
raw_data$Ind_ID <- as.factor(raw_data$Ind_ID)

# Make a raw data sheet with no missing length & mass data
raw_data_lengthmass <- filter(raw_data, !is.na(Mass) & !is.na(Length))
raw_data_lengthmass$treatment_combined <- droplevels(raw_data_lengthmass$treatment_combined)
raw_data_lengthmass <- raw_data_lengthmass %>% 
  mutate(Mass_c = Mass - mean(Mass),
         Length_c = Length - mean(Length))

# Take a look at the raw data
ggplot(data = raw_data, aes(x = Time, y = ATPase)) +
  geom_point() +
  facet_wrap(~Treatment)

ggplot(data = raw_data %>% group_by(Ind_ID) %>% slice(1), aes(x = rep(1), y = Length)) +
  geom_violin() +
  geom_point(size = 2) +
  facet_wrap(~Treatment) +
  theme_bw(base_size = 16)

ggplot(data = raw_data %>% group_by(Ind_ID) %>% slice(1), aes(x = rep(1), y = Mass)) +
  geom_violin() +
  geom_point(size = 2) +
  facet_wrap(~Treatment) +
  theme_bw(base_size = 16)

ggplot(data = raw_data %>% group_by(Ind_ID) %>% slice(1), aes(x = Time, y = Length)) +
  geom_violin() +
  geom_point(size = 2) +
  facet_wrap(~Treatment) +
  theme_bw(base_size = 16)

ggplot(data = raw_data %>% group_by(Ind_ID) %>% slice(1), aes(x = Time, y = Mass)) +
  geom_violin() +
  geom_point(size = 2) +
  facet_wrap(~Treatment) +
  theme_bw(base_size = 16)

# Run the model. -1 removed the intercept, Run is used as a random effect, and the treatment and timepoint variables were combined here
# Use family = skew_normal because when using family = gaussian, the posterior distribution was systematically left-skewed
# Create three models to test how much length and mass information contributes to model fit. 
# First, mass interacting with length
lengthmass_model <- brm(data = raw_data_lengthmass, family = skew_normal,
                 ATPase ~ -1 + (1 | Run) + (1 | Ind_ID) + treatment_combined + (Mass_c * Length_c), 
                 prior = c(prior(uniform(0, 80), lb = 0, ub = 80, class = b),
                           prior(cauchy(4, 5), class = sigma)),
                 iter = 5000, warmup = 1000, chains = 4, cores = 4,  
                 control = list(adapt_delta = 0.999, max_treedepth = 20))

# Length only with no mass
length_model <- brm(data = raw_data_lengthmass, family = skew_normal,
                        ATPase ~ -1 + (1 | Run) + (1 | Ind_ID) + treatment_combined + Length_c, 
                        prior = c(prior(uniform(0, 80), lb = 0, ub = 80, class = b),
                                  prior(cauchy(4, 5), class = sigma)),
                        iter = 5000, warmup = 1000, chains = 4, cores = 8,  
                        control = list(adapt_delta = 0.98, max_treedepth = 20))

# No mass or length
nomorphology_model <- brm(data = raw_data_lengthmass, family = skew_normal,
                    ATPase ~ -1 + (1 | Run) + (1 | Ind_ID) + treatment_combined, 
                    prior = c(prior(uniform(0, 80), lb = 0, ub = 80, class = b),
                              prior(cauchy(4, 5), class = sigma)),
                    iter = 5000, warmup = 1000, chains = 4, cores = 8,  
                    control = list(adapt_delta = 0.98, max_treedepth = 20))

# Investigate WAIC by adding it as a criterion to the different models
waic_lengthmass <- add_criterion(lengthmass_model, "waic")
waic_length <- add_criterion(length_model, "waic")
waic_nomorphology <- add_criterion(nomorphology_model, "waic")

waic_diff <- loo_compare(waic_lengthmass, waic_length, waic_nomorphology, criterion = "waic") 
print(waic_diff, simplify = FALSE)

# Use BayesTestR to look at Bayes factors among the models
bayesfactor_models(lengthmass_model, length_model, nomorphology_model, denominator = 1)

# Plot WAIC values
waic_diff[, 7:8] %>% 
  data.frame() %>% 
  rownames_to_column(var = "model_name") %>% 
  
  ggplot(aes(x    = model_name, 
             y    = waic, 
             ymin = waic - se_waic, 
             ymax = waic + se_waic)) +
  geom_pointrange(shape = 21, color = carto_pal(7, "BurgYl")[7], fill = carto_pal(7, "BurgYl")[5]) +
  coord_flip() +
  labs(y = "Watanabe-Akaike Information Criterion", x = NULL,
       title = "Comparing Different Models") +
  theme_classic() +
  theme(text             = element_text(family = "Courier"),
        axis.ticks.y     = element_blank(),
        panel.background = element_rect(fill = alpha(carto_pal(7, "BurgYl")[3], 1/4)))

# WAIC was lowest for the model with length and mass interacting, but just barely. Therefore, here is a model with the full dataset including the 720h fish without length and mass. While those variables cannot be included, model fit was marginally affected and the extra timepoint was informative
fulldata_model <- brm(data = raw_data, family = skew_normal,
                        bf(ATPase ~ treatment_combined + (1 | Run) + (1 | Ind_ID),
                           sigma ~ treatment_combined), 
                        prior = c(prior(normal(0, 20), class = Intercept),
                                  prior(normal(0, 20), class = b),
                                  prior(cauchy(0,4), class = sd)),
                        iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                        control = list(adapt_delta = 0.95, max_treedepth = 14))
# prior(cauchy(4, 5), class = sigma)
#pairs(brm_model)
# Take a look at the overall model
summary(fulldata_model)

# Use BayesTestR to look at posterior and priors
describe_posterior(fulldata_model)
describe_prior(fulldata_model)
check_prior(fulldata_model)
p_map(fulldata_model)

# Look at trace plots
plot(fulldata_model, N = 3, ask = FALSE)

# Look at the posterior distribution. It was left skewed when using a Gaussian family distribution, but was fixed by using skew_normal instead
pp_check(fulldata_model, ndraws = 100) +
  theme_bw(base_size = 20)


# Plot conditional effects of treatment on ATPase activity
plot(conditional_effects(fulldata_model), points = T)

# Use gather emmeans draws to plot the model results, showing different posterior distributions for each treatment at each timepoint

atpase_emmeansdraws <- gather_emmeans_draws(emmeans(fulldata_model, ~ treatment_combined)) |>
  separate(treatment_combined, into = c("treatment", "time"), sep = "_") |> 
  mutate(treatment = factor(treatment, levels = c("Combo", "pCO2", "Temp", "Control"))) |> 
  mutate(time = factor(time, levels = c("0", "6", "24", "168", "720"))) 


atpase_plot <- ggplot(data = atpase_emmeansdraws, aes(x = .value, y = time, fill = treatment, colour = treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.5) +
  #geom_point(data = raw_data, aes(x = ATPase, y = Time, fill = Treatment, colour = Treatment), alpha = 0.5, shape = 15) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  xlab(bquote(Na^"+"*"/"*K^"+"~ATPase~Activity~(nmol~NADH*"\U00B7"*min^-1*"\U00B7"*mg~pr^-1))) +
  ylab("Timepoint (hours)") +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment")) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.91, 0.14), legend.background = element_blank())
atpase_plot
ggsave(filename = "atpase_plot.pdf", plot = atpase_plot, dpi = 2000, width = 14, height = 8.5)

save(atpase_plot, atpase_emmeansdraws, file = "atpase_plot.RData")



# Pull pairwise results from the emmeans function with tidybayes, and make halfeye plots
atpase_pairwise_table <- pairs(emmeans(fulldata_model, ~ treatment_combined))
atpase_pairwise_table <- summary(atpase_pairwise_table)

write_tsv(x = atpase_pairwise_table, file = "atpase_pairwise_contrasts.txt")

pairwise_ATPase_plot <- gather_emmeans_draws(pairs(emmeans(fulldata_model, ~ treatment_combined))) %>% 
  ggplot(aes(x = .value, y = contrast)) + 
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab(bquote(Estimated~Pairwise~Difference~"in"~Na^"+"*"/"*K^"+"~ATPase~Activity~(nmol~NADH*"\U00B7"*min^-1*"\U00B7"*mg*"\U00B7"*pr^-1))) +
  ylab("Treatment & Timepoint") +
  theme_bw()

ggsave(filename = "pairwise_atpase_plot.pdf", plot = pairwise_ATPase_plot, dpi = 2000, width = 14, height = 20)

# Remove the pairwise plot since it is a huge gg object
rm(pairwise_ATPase_plot)
#########################################################################
# The code below was used for fun but was not necessary for publication

# Test different hypotheses with the data
plot(hypothesis(fulldata_model, "treatment_combinedpCO2_0 > treatment_combinedControl_0", alpha = 0.05))
hypothesis(fulldata_model, "treatment_combinedpCO2_0 > treatment_combinedControl_0", alpha = 0.05)

hypothesis(fulldata_model, "treatment_combinedpCO2_0 + treatment_combinedpCO2_24 + treatment_combinedpCO2_168 + treatment_combinedpCO2_720 > treatment_combinedControl_0 + treatment_combinedControl_24 + treatment_combinedControl_168 + treatment_combinedControl_720", alpha = 0.05)

plot(hypothesis(fulldata_model, "treatment_combinedpCO2_0 + treatment_combinedpCO2_24 + treatment_combinedpCO2_168 + treatment_combinedpCO2_720 > treatment_combinedControl_0 + treatment_combinedControl_24 + treatment_combinedControl_168 + treatment_combinedControl_720", alpha = 0.05))

hypothesis(fulldata_model, "treatment_combinedTemp_168 > treatment_combinedTemp_720", alpha = 0.05)

#############################################
# Use DHARMa to check how well the model was fit
# sample from the Posterior Predictive Distribution
# Basic model residual checks using the notes at https://github.com/florianhartig/DHARMa/issues/33
preds <- t(posterior_predict(fulldata_model))

res <- createDHARMa(simulatedResponse = preds, 
                    fittedPredictedResponse = apply(preds, 1, mean), 
                    observedResponse = raw_data$ATPase, 
                    integerResponse = F)

plot(res, quantreg = FALSE)
rm(preds, res)

# A purpose-built function by F. Rodriguez-Sanchez for checking brms models with DHARMa. Following his guide:
# https://frodriguezsanchez.net/post/using-dharma-to-check-bayesian-models-fitted-with-brms/
check_brms <- function(model,             # brms model
                       integer = FALSE,   # integer response? (TRUE/FALSE)
                       plot = TRUE,       # make plot?
                       ...                # further arguments for DHARMa::plotResiduals 
) {
  
  mdata <- brms::standata(model)
  if (!"Y" %in% names(mdata))
    stop("Cannot extract the required information from this brms model")
  
  dharma.obj <- DHARMa::createDHARMa(
    simulatedResponse = t(brms::posterior_predict(model, ndraws = NULL)),
    observedResponse = mdata$Y, 
    fittedPredictedResponse = apply(
      t(brms::posterior_epred(model, nsamples = 1000, re.form = NA)),
      1,
      mean),
    integerResponse = integer)
  
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  
  invisible(dharma.obj)
  
}

# Run the model check function, then test for dispersion
model_check <- check_brms(brm_model, integer = FALSE)
testDispersion(model_check)

#####################################




# Look at a plot of predicted ATPase activity over time
time_plot <- raw_data %>% 
  add_predicted_draws(fulldata_model) %>% 
  ggplot(aes(x = Time, y = ATPase))  +
  stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), color = "#08519C") +
  geom_point(data = raw_data) +
  scale_fill_brewer() +
  facet_wrap(~ Treatment) +
  theme_bw()
#time_plot
ggsave(filename = "time_plot.pdf", plot = time_plot, dpi = 900)
rm(time_plot) # Remove the time plot after exporting it because it's a huge R object

# Look at the variables available in the model
get_variables(fulldata_model)

# Pull variables. Here, treatments all start with "b" so use regex to pull just those
model_dataframe <- fulldata_model %>% 
  spread_draws(., `[\\(b].*`, regex = TRUE) %>% 
  pivot_longer(cols = b_treatment_combinedCombo_0:b_treatment_combinedTemp_720, names_to = "treatment", values_to = "est") %>% 
  separate(treatment, into = c("b", "treat", "treatment", "time"), sep = "_") %>% 
  select(-b, -treat) %>% 
  mutate(treatment = str_remove(treatment, "combined"))
model_dataframe$treatment <- factor(model_dataframe$treatment, levels = c("Combo", "pCO2", "Temp", "Control"))
model_dataframe$time <- factor(model_dataframe$time, levels = c("0", "6", "24", "168", "720"))

# Create a smaller version of the model dataframe just for troubleshooting plots
#subset_results <- model_dataframe %>% 
  #group_by(treatment, time) %>% 
  #slice_sample(n = 500)

# Make halfeye plots of the different treatments
atpase_plot <- ggplot(data = model_dataframe, aes(x = est, y = time, fill = treatment, colour = treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.5) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  xlab(bquote(Na^"+"*"/"*K^"+"~ATPase~Activity~(nmol~NADH*"\U00B7"*min^-1*"\U00B7"*mg*"\U00B7"*pr^-1))) +
  ylab("Timepoint (hours)") +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment")) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.91, 0.14), legend.background = element_blank())
atpase_plot

# Save the halfeye plot
ggsave(filename = "atpase_plot.pdf", plot = atpase_plot, dpi = 2000, width = 14, height = 8.5)

# a second version of the plot with less emphasis on the distribution tails
atpase_plot_gradient <- ggplot(data = model_dataframe, aes(x = est, y = time, fill = treatment, colour = treatment)) +
  stat_slab(aes(alpha = 0.75, fill = treatment, fill_ramp = stat(f), slab_linetype = NA)) +
  ggdist::scale_fill_ramp_continuous() +
  stat_pointinterval(alpha = 0.75, show.legend = FALSE) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  xlab(bquote(Na^"+"*"/"*K^"+"~ATPase~Activity~(nmol~NADH*"\U00B7"*min^-1*"\U00B7"*mg*"\U00B7"*pr^-1))) +
  ylab("Timepoint (hours)") +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment"), alpha = "none", fill_ramp = "none") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.91, 0.14), legend.background = element_blank())

#ggsave(filename = "atpase_plot_gradient.pdf", plot = atpase_plot_gradient, dpi = 2000, width = 14, height = 8.5)
ggsave(filename = "atpase_plot_gradient.png", plot = atpase_plot_gradient, dpi = 2000, width = 14, height = 8.5)


############################################
# Use performance to do other model checks

# Look at different credible intervals and among the different treatments
posterior_intervals <- describe_posterior(brm_model, test = c("pd", "ROPE", "BF"))

# Marginal R2 provides the variance explained only by fixed effects and conditional R2 provides the variance explained by the entire model, i.e., both fixed effects and random effects.
model_performance <- model_performance(brm_model)



#######################################
# Use modelbased to visualize the model, following this guide
# https://easystats.github.io/bayestestR/articles/example2.html
vizdata <- estimate_relation(fulldata_model)

ggplot(vizdata, aes(x = treatment_combined, y = Predicted)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), alpha = 0.5) +
  ylab("Predicted ATPase Activity") +
  theme_bw()

#######################################
# Make a plot of all pairwise differences
pairwise_plot <- fulldata_model %>% 
  spread_draws(., `[\\(b].*`, regex = TRUE) %>% 
  pivot_longer(cols = b_treatment_combinedCombo_0:b_treatment_combinedTemp_720, names_to = "treatment", values_to = "est") %>% 
  mutate(treatment = as.factor(treatment), est = as.numeric(est)) %>% 
  mutate(treatment = str_remove(treatment, "b_treatment_combined")) %>% 
  compare_levels(est, by = treatment) %>% 
  ggplot(aes(y = treatment, x = est)) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Estimated Pairwise Difference") +
  ylab("Treatment") +
  theme_bw()

ggsave(filename = "pairwise_ATPase.pdf", plot = pairwise_plot, dpi = 900, height = 24, width = 8)


# Pull particular differences that we care about and plot those
select_pairwise_draws <- fulldata_model %>% 
  spread_draws(., `[\\(b].*`, regex = TRUE) %>% 
  pivot_longer(cols = b_treatment_combinedCombo_0:b_treatment_combinedTemp_720, names_to = "treatment", values_to = "est") %>% 
  mutate(treatment = as.factor(treatment), est = as.numeric(est)) %>% 
  mutate(treatment = str_remove(treatment, "b_treatment_combined")) %>% 
  compare_levels(est, by = treatment, comparison = list(c("Control_0", "Control_6"),
                                                        c("Control_6", "Control_24"),
                                                        c("Control_24", "Control_168"),
                                                        c("Control_168", "Control_720"),
                                                        c("Temp_0", "Temp_6"),
                                                        c("Temp_6", "Temp_24"),
                                                        c("Temp_24", "Temp_168"),
                                                        c("Temp_168", "Temp_720"))) %>% 
  mutate(treatment = factor(treatment, levels = c("Control_0 - Control_6", 
                                                  "Control_6 - Control_24",
                                                  "Control_24 - Control_168",
                                                  "Control_168 - Control_720",
                                                  "Temp_0 - Temp_6",
                                                  "Temp_6 - Temp_24",
                                                  "Temp_24 - Temp_168",
                                                  "Temp_168 - Temp_720"))) 

select_pairwise_plot <- ggplot(select_pairwise_draws, aes(y = treatment, x = est, group = treatment)) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Estimated Pairwise Difference") +
  ylab("Treatment") +
  theme_bw()

ggsave(filename = "select_pairwise_ATPase.pdf", plot = select_pairwise_plot, dpi = 900, height = 8.5, width = 12, scale = 0.5)


