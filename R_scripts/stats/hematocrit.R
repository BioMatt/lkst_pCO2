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
library(emmeans) # Estimated marginal means for pulling draws and pairwise comparisons
# Faster Hamiltonian Monte Carlo sampling
options(brms.backend = "cmdstanr")

raw_data <- readxl::read_excel("../sturgeon blood work May 2018_MT.xlsx", sheet = "hematocrit") %>% 
  separate(Ind, into = c("Time", "Index")) %>% 
  mutate(Time = case_when(Time == "Recovery" ~ "month",
                          Time == "0d" ~ "zero", 
                          Time == "24h" ~ "day",
                          Time == "168h" ~ "week",
                          Time == "6h" ~ "six_hours")) %>% 
  mutate(Treatment = case_when(Treatment == "pCO2+Temp" ~ "pCO2.Temp",
                               Treatment == "Control" ~ "Control",
                               Treatment == "pCO2" ~ "pCO2",
                               Treatment == "Temp" ~ "Temp")) %>% 
  mutate(treatment_combined = paste0(Treatment, "_", Time))  %>% 
  mutate(treatment_combined = factor(treatment_combined))

# Take a look at the raw data
ggplot(data = raw_data, aes(x = Time, y = Hematocrit)) +
  geom_point() +
  facet_wrap(~Treatment)

# Run the model
hematocrit_model <- brm(data = raw_data, family = skew_normal,
                      bf(Hematocrit ~ treatment_combined,
                         sigma ~ treatment_combined), 
                      prior = c(prior(normal(10, 30), class = Intercept),
                              prior(normal(10, 30), class = b)),
                      iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                      control = list(adapt_delta = 0.99, max_treedepth = 12))

# Take a look at the overall model
summary(hematocrit_model)

# Look at trace plots
plot(hematocrit_model, N = 2, ask = FALSE)

# Look at the posterior distribution. It was left skewed when using a Gaussian family distribution, but was fixed by using skew_normal instead
pp_check(hematocrit_model, ndraws = 100) +
  theme_bw(base_size = 20)

plot(conditional_effects(hematocrit_model), points = TRUE)

##################################################
# Use gather emmeans draws to plot the model results, showing different posterior distributions for each treatment at each timepoint

emmeansdraws <- gather_emmeans_draws(emmeans(hematocrit_model, ~ treatment_combined, re_formula = NULL)) |>
  separate(treatment_combined, into = c("treatment", "time"), sep = "_") |> 
  mutate(treatment = factor(treatment, levels = c("pCO2.Temp", "pCO2", "Temp", "Control"))) |> 
  mutate(Time = factor(time, levels = c("zero", "six_hours", "day", "week", "month"))) 


hematocrit_plot <- ggplot(data = emmeansdraws, aes(x = .value, y = time, fill = treatment, colour = treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.5) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  xlab("Hematocrit (%)") +
  ylab("Timepoint (hours)") +
  scale_y_discrete(labels = c("0", "6", "24", "168", "720")) +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment")) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.91, 0.14), legend.background = element_blank())
#hematocrit_plot

ggsave(filename = "hematocrit_plot.pdf", plot = hematocrit_plot, dpi = 2000, width = 14, height = 8.5)



# Pull pairwise results from the emmeans function with tidybayes, and make halfeye plots
pairwise_table <- pairs(emmeans(hematocrit_model, ~ treatment_combined, re_formula = NULL))
pairwise_table <- summary(pairwise_table)

write_tsv(x = pairwise_table, file = "hematocrit_pairwise_contrasts.txt")

pairwise_hematocrit_plot <- gather_emmeans_draws(pairs(emmeans(hematocrit_model, ~ treatment_combined, re_formula = NULL))) %>% 
  ggplot(aes(x = .value, y = contrast)) + 
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Estimated Pairwise Difference in Hematocrit") +
  ylab("Treatment & Timepoint") +
  theme_bw()

ggsave(filename = "pairwise_hematocrit_plot.pdf", plot = pairwise_hematocrit_plot, dpi = 2000, width = 14, height = 20)

# Remove the pairwise plot since it is a huge gg object
rm(pairwise_hematocrit_plot)

#############################################

# Use BayesTestR to look at posterior and priors
describe_posterior(hematocrit_model)
describe_prior(hematocrit_model)
check_prior(hematocrit_model)
p_map(hematocrit_model)

# Plot conditional effects of treatment on hematocrit
plot(conditional_effects(hematocrit_model), points = T)

# Look at the variables available in the model
get_variables(hematocrit_model)

# Pull variables. Here, treatments all start with "b" so use regex to pull just those
model_dataframe <- hematocrit_model %>% 
  spread_draws(., `[\\(b].*`, regex = TRUE) %>% 
  pivot_longer(cols = b_treatment_combinedControl_168:b_treatment_combinedTemp_720, names_to = "treatment", values_to = "est") %>% 
  separate(treatment, into = c("b", "treat", "treatment", "time"), sep = "_") %>% 
  select(-b, -treat) %>% 
  mutate(treatment = str_remove(treatment, "combined"))
model_dataframe$treatment <- factor(model_dataframe$treatment, levels = c("pCO2.Temp", "pCO2", "Temp", "Control"))
model_dataframe$time <- factor(model_dataframe$time, levels = c("0", "6", "24", "168", "720"))


hematocrit_plot <- ggplot(data = model_dataframe, aes(x = est, y = time, fill = treatment, colour = treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.5) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  xlab("Hematocrit (% RBC)") +
  ylab("Timepoint (hours)") +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment")) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.91, 0.14), legend.background = element_blank())
hematocrit_plot

# Save the halfeye plot
ggsave(filename = "hematocrit_plot.pdf", plot = hematocrit_plot, dpi = 2000, width = 14, height = 8)

