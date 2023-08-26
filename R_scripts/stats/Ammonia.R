library(tidyverse) # Handling large datasets
library(brms) # Run the Bayesian models
library(tidybayes) # Plot Bayesian models
library(modelr) # Make model functions work in pipes
library(fishualize) # Get nice fish-based colours
library(DHARMa) # Check residuals and dispersion for Bayesian models
library(performance) # Check model performance
library(modelbased) # Visualize the model's posterior predictions
library(emmeans) # Estimated marginal means for pulling draws and pairwise comparisons
# Faster Hamiltonian Monte Carlo sampling
options(brms.backend = "cmdstanr")

# Read in the raw data. Separate fish IDs from timepoint information, replace timepoint information with the hours represented by it. Pivot the dataframe longer for modeling, take out + signs in the treatment names, create a new column of combined treatment and timepoint, turn the combined column into a factor, then re-order the table.
raw_data <- readxl::read_excel("../Copy of sturgeon MO2 and ammonia excretion 2018_MT.xlsx", sheet = "Ammonia") %>% 
  separate(Fish_ID, into = c("Timepoint", "ID")) %>% 
  mutate(Timepoint = case_when(Timepoint == "0d" ~ 0,
                               Timepoint == "24h" ~ 24,
                               Timepoint == "168h" ~ 168,
                               Timepoint == "Recovery" ~ 720)) %>% 
  pivot_longer(cols = Control:`pCO2+Temp`, names_to = "Treatment", values_to = "Ammonia") %>%
  mutate(Treatment = gsub("\\+", ".", Treatment)) %>% 
  mutate(Treatment_Combined = paste0(Treatment, "_", Timepoint)) %>% 
  mutate(Treatment_Combined = as.factor(Treatment_Combined)) %>% 
  select(ID, Timepoint, Treatment, Treatment_Combined, Ammonia) %>% 
  filter(!is.na(Ammonia))

# Take a look at the raw data
ggplot(data = raw_data, aes(x = Timepoint, y = Ammonia)) +
  geom_point() +
  facet_wrap(~Treatment)

# Run the model. -1 removed the intercept, Run is used as a random effect, and the treatment and timepoint variables were combined here
# Use family = skew_normal because when using family = gaussian, the posterior distribution was systematically left-skewed
brm_model <- brm(data = raw_data, family = skew_normal,
                 bf(Ammonia ~ Treatment_Combined,
                    sigma ~ Treatment_Combined), 
                 prior = c(prior(normal(0, 1), class = b),
                           prior(normal(0, 1), class = Intercept)),
                 iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                 control = list(adapt_delta = 0.98))


#pairs(brm_model)
# Take a look at the overall model
summary(brm_model)

# Look at trace plots
plot(brm_model, N = 3, ask = FALSE)

# Look at the posterior distribution. It was left skewed when using a Gaussian family distribution, but was fixed by using skew_normal instead
pp_check(brm_model, ndraws = 100) +
  theme_bw(base_size = 20)

# Look at conditional effects
plot(conditional_effects(brm_model), points = TRUE)

##################################################
# Use gather emmeans draws to plot the model results, showing different posterior distributions for each treatment at each timepoint

emmeansdraws <- gather_emmeans_draws(emmeans(brm_model, ~ Treatment_Combined)) |>
  separate(Treatment_Combined, into = c("treatment", "time"), sep = "_") |> 
  mutate(treatment = factor(treatment, levels = c("pCO2.Temp", "pCO2", "Temp", "Control"))) |> 
  mutate(time = factor(time, levels = c("0", "24", "168", "720"))) 


ammonia_plot <- ggplot(data = emmeansdraws, aes(x = .value, y = time, fill = treatment, colour = treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.5) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  xlab(bquote(Ammonia~Excretion~(mu*mol*"\U00B7"*g^-1*"\U00B7"*h^-1))) +
  ylab("Timepoint (hours)") +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment")) +
  xlim(0.0, 1.1) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.91, 0.14), legend.background = element_blank())

ggsave(filename = "ammonia_plot.pdf", plot = ammonia_plot, dpi = 2000, width = 14, height = 8.5)



# Pull pairwise results from the emmeans function with tidybayes, and make halfeye plots
pairwise_table <- pairs(emmeans(brm_model, ~ Treatment_Combined))
pairwise_table <- summary(pairwise_table)

write_tsv(x = pairwise_table, file = "ammonia_pairwise_contrasts.txt")

pairwise_ammonia_plot <- gather_emmeans_draws(pairs(emmeans(brm_model, ~ Treatment_Combined))) %>% 
  ggplot(aes(x = .value, y = contrast)) + 
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab(bquote(Estimated~Pairwise~Difference~"in"~Ammonia~Excretion~(mu*mol*"\U00B7"*g^-1*"\U00B7"*h^-1))) +
  ylab("Treatment & Timepoint") +
  theme_bw()

ggsave(filename = "pairwise_ammonia_plot.pdf", plot = pairwise_ammonia_plot, dpi = 2000, width = 14, height = 20)

# Remove the pairwise plot since it is a huge gg object
rm(pairwise_ammonia_plot)

#############################################
# Use DHARMa to check how well the model was fit
# sample from the Posterior Predictive Distribution
# Basic model residual checks using the notes at https://github.com/florianhartig/DHARMa/issues/33
preds <- t(posterior_predict(brm_model))

res <- createDHARMa(simulatedResponse = preds, 
                    fittedPredictedResponse = apply(preds, 1, mean), 
                    observedResponse = raw_data$Ammonia, 
                    integerResponse = F)

plot(res, quantreg = FALSE)

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
# Plot conditional effects of treatment on Ammonia activity
plot(conditional_effects(brm_model), points = T)

# Look at the variables available in the model
get_variables(brm_model)

# Pull variables. Here, treatments all start with "b" so use regex to pull just those
model_dataframe <- brm_model %>% 
  spread_draws(., `[\\(b].*`, regex = TRUE) %>% 
  pivot_longer(cols = b_Treatment_CombinedControl_0:b_Treatment_CombinedTemp_720, names_to = "treatment", values_to = "est") %>% 
  separate(treatment, into = c("b", "treat", "treatment", "time"), sep = "_") %>% 
  select(-b, -treat) %>% 
  mutate(treatment = str_remove(treatment, "Combined"))
model_dataframe$treatment <- factor(model_dataframe$treatment, levels = c("pCO2.Temp", "pCO2", "Temp", "Control"))
model_dataframe$time <- factor(model_dataframe$time, levels = c("0", "6", "24", "168", "720"))

# Create a smaller version of the model dataframe just for troubleshooting plots
subset_results <- model_dataframe %>% 
  group_by(treatment, time) %>% 
  slice_sample(n = 500)


ammonia_plot <- ggplot(data = model_dataframe, aes(x = est, y = time, fill = treatment, colour = treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.5) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  xlab(bquote(Ammonia~Excretion~(mu*mol*"\U00B7"*g^-1*"\U00B7"*h^-1))) +
  ylab("Timepoint (hours)") +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment")) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.91, 0.14), legend.background = element_blank())

ggsave(filename = "ammonia_plot.pdf", plot = ammonia_plot, dpi = 2000, width = 14, height = 8.5)
