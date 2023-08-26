library(tidyverse) # Handling large datasets
library(brms) # Run the Bayesian models
library(tidybayes) # Plot Bayesian models
library(modelr) # Make model functions work in pipes
library(fishualize) # Get nice fish-based colours
library(performance) # Check model performance
library(modelbased) # Visualize the model's posterior predictions
library(emmeans) # Estimated marginal means to compare groups

# Faster Hamiltonian Monte Carlo sampling
options(brms.backend = "cmdstanr")


# Read in the raw data
raw_data <- readxl::read_excel("2022 Luke MO2 30 day recovery.xlsx", sheet = 1) %>% 
  mutate(Date = factor(Date)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("control", "temp", "pco2", "combo")))


# Run a model of SMR but controlled for umol/g/h

# Take a look at the raw data
ggplot(data = raw_data, aes(x = Treatment, y = RMR)) +
  geom_point()
ggplot(data = raw_data, aes(x = Mass*Length, y = RMR)) +
  geom_point()

brm_RMR <- brm(data = raw_data, family = skew_normal,
                 bf(RMR ~ Treatment + (1 | Date),
                    sigma ~ Treatment), 
                 prior = c(prior(normal(50, 50), class = b),
                           prior(normal(50, 50), class = Intercept)),
                 iter = 20000, warmup = 5000, chains = 4, cores = 4,
                 control = list(adapt_delta = 0.98, max_treedepth = 14))


#pairs(brm_model)
# Take a look at the overall model
summary(brm_RMR)

# Look at trace plots
plot(brm_RMR, N = 3, ask = FALSE)

# Look at the posterior distribution
pp_check(brm_RMR, ndraws = 100) +
  theme_bw(base_size = 20)

# Look at the variables available in the model
get_variables(brm_RMR)

plot(conditional_effects(brm_RMR), points = TRUE, theme = theme_b())
# Pull variables. Here, treatments all start with "b" so use regex to pull just those
# Code commented out because gather emmeans draws turned out to be a better way to get model posteriors
# model_dataframe_RMR <- brm_RMR %>% 
#   spread_draws(., `[\\(b].*`, regex = TRUE) %>% 
#   pivot_longer(cols = b_Treatmentcombo:b_Treatmenttemp, names_to = "treatment", values_to = "est") %>% 
#   separate(treatment, into = c("b", "treatment"), sep = "_") %>% 
#   select(-b) %>% 
#   mutate(treatment = str_remove(treatment, "Treatment")) %>% 
#   mutate(treatment = factor(treatment, levels = c("control", "temp", "pco2", "combo")))

rmr_draws <- gather_emmeans_draws(emmeans(brm_RMR, ~ Treatment)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("control", "temp", "pco2", "combo")))


RMR_plot <- ggplot(data = rmr_draws, aes(x = .value, y = Treatment, fill = Treatment, colour = Treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.8, show.legend = FALSE) +
  labs(x = bquote('Routine Metabolic Rate ('*mg~O[2]*"\U00B7"*kg^-1*"\U00B7"*Hr^-1*')'), y = "Treatment") +
  scale_y_discrete(labels = c("Control", "Temperature", bquote(italic(p)*CO[2]), bquote(italic(p)*CO[2]*"+"*Temp))) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = 1, labels = rev(c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control"))) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = 1, labels = rev(c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control"))) +
  guides(fill = guide_legend(title="Treatment", reverse = TRUE), colour = guide_legend(title="Treatment", reverse = TRUE)) +
  theme_bw(base_size = 18)

RMR_plot

ggsave(filename = "RMR_plot.pdf", plot = RMR_plot, dpi = 2000, scale = 1.5)

###################################################
# The same type of model but for MMR
# Take a look at the raw data
ggplot(data = raw_data, aes(x = Treatment, y = MMR)) +
  geom_point()


brm_MMR <- brm(data = raw_data, family = skew_normal,
               bf(MMR ~ Treatment + (1 | Date),
                  sigma ~ Treatment), 
               prior = c(prior(normal(50, 50), class = b),
                         prior(normal(50, 50), class = Intercept)),
               iter = 50000, warmup = 10000, chains = 4, cores = 4,
               control = list(adapt_delta = 0.99, max_treedepth = 14))

#pairs(brm_model)
# Take a look at the overall model
summary(brm_MMR)

# Look at trace plots
plot(brm_MMR, N = 3, ask = FALSE)

# Look at the posterior distribution
pp_check(brm_MMR, ndraws = 100) +
  theme_bw(base_size = 20)

# Look at the variables available in the model
get_variables(brm_MMR)

# Pull variables. Here, treatments all start with "b" so use regex to pull just those
# model_dataframe_MMR <- brm_MMR %>% 
#   spread_draws(., `[\\(b].*`, regex = TRUE) %>% 
#   pivot_longer(cols = b_Treatmentcombo:b_Treatmenttemp, names_to = "treatment", values_to = "est") %>% 
#   separate(treatment, into = c("b", "treatment"), sep = "_") %>% 
#   select(-b) %>% 
#   mutate(treatment = str_remove(treatment, "Treatment")) %>% 
#   mutate(treatment = factor(treatment, levels = c("control", "temp", "pco2", "combo")))

mmr_draws <- gather_emmeans_draws(emmeans(brm_MMR, ~ Treatment)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("control", "temp", "pco2", "combo")), Time = rep("720h", 1))


MMR_plot <- ggplot(data = mmr_draws, aes(x = .value, y = Treatment, fill = Treatment, colour = Treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.8, show.legend = FALSE) +
  labs(x = bquote('Maximum Metabolic Rate ('*mg~O[2]*"\U00B7"*kg^-1*"\U00B7"*Hr^-1*')'), y = "Treatment") +
  scale_y_discrete(labels = c("Control", "Temperature", bquote(italic(p)*CO[2]), bquote(italic(p)*CO[2]*"+"*Temp))) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = 1, labels = rev(c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control"))) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = 1, labels = rev(c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control"))) +
  guides(fill = guide_legend(title="Treatment", reverse = TRUE), colour = guide_legend(title="Treatment", reverse = TRUE)) +
  theme_bw(base_size = 18)
MMR_plot

ggsave(filename = "MMR_plot.pdf", plot = MMR_plot, dpi = 2000, scale = 1.5)

MMR_plot2 <- ggplot(data = mmr_draws, aes(x = .value, y = Time, fill = Treatment, colour = Treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.5) +
  labs(x = bquote('Maximum Metabolic Rate ('*mg~O[2]*"\U00B7"*kg^-1*"\U00B7"*Hr^-1*')'), y = "Treatment") +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = rev(c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control"))) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = rev(c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control"))) +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment")) +
  theme_bw(base_size = 18)
MMR_plot2
# Take a look at the pCO2 versus combined fish in the MMR model
hypothesis(brm_MMR, "Treatmentpco2 > Treatmentcombo", alpha = 0.05)
plot(hypothesis(brm_MMR, "Treatmentpco2 > Treatmentcombo", alpha = 0.05))

######################################################################################
# Use emmeans as a similar method for comparing groups in a pairwise way. This ended up being simpler because each contrast did not need to be specified
# Start with RMR
# Look at the means between populations
emmeans(brm_RMR, ~ Treatment)
pairs(emmeans(brm_RMR, ~ Treatment))
write_tsv(file = "RMR_emmeans.txt", x = summary(pairs(emmeans(brm_RMR, ~ Treatment))))

# Plot the emmeans credible intervals without distributions
plot(pairs(emmeans(brm_RMR, ~ Treatment))) +
  theme_bw()

# Pull results from the emmeans pairs function with tidybayes, and make halfeye plots
RMR_pairwise_plot_emmeans <- gather_emmeans_draws(pairs(emmeans(brm_RMR, ~ Treatment))) |>
  ggplot(aes(x = .value, y = contrast)) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Estimated Pairwise Difference in Routine Metabolic Rate") +
  ylab("Contrast") +
  theme_bw(base_size = 8)

ggsave(filename = "RMR_pairwise.pdf", plot = RMR_pairwise_plot_emmeans, dpi = 900)

# Also check the specific comparison of pCO2+temp against pCO2 with hypothesis
hypothesis(x = brm_RMR, hypothesis = "Treatmentpco2 > Treatmentcombo")
hypothesis(x = brm_RMR, hypothesis = "Treatmenttemp > Treatmentcombo")


# Now MMR
# Look at the means between populations
emmeans(brm_MMR, ~ Treatment)
write_tsv(file = "MMR_emmeans.txt", x = summary(pairs(emmeans(brm_MMR, ~ Treatment))))

# Plot the emmeans credible intervals without distributions
plot(pairs(emmeans(brm_MMR, ~ Treatment))) +
  theme_bw()

# Pull results from the emmeans pairs function with tidybayes, and make halfeye plots
MMR_pairwise_plot_emmeans <- gather_emmeans_draws(pairs(emmeans(brm_MMR, ~ Treatment))) |>
  ggplot(aes(x = .value, y = contrast)) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Estimated Pairwise Difference in Maximum Metabolic Rate") +
  ylab("Contrast") +
  theme_bw(base_size = 8)

ggsave(filename = "MMR_pairwise.pdf", plot = MMR_pairwise_plot_emmeans, dpi = 900)


# An unrelated plot for visualizing normal distributions and credible intervals
test <- ggplot(data = model_dataframe_RMR, aes(x = est, y = treatment)) +
  stat_halfeye(orientation = "horizontal", fill = "#2b8cbe") +
  theme_void()

ggsave(filename = "test.pdf", plot = test, dpi = 3000)
