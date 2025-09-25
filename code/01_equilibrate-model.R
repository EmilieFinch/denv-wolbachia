#### Code to calibrate serotype or strain-level compartmental  models ####
## Equilibrate the model for 100 years using static demography using demographic data for Brazil from 1974 until 2024 
## and assuming static demography prior to this point

source(here("code", "utils.R"))

eq_type <- "strain" # can be strain or serotype

### Set paths
output_path <- here("output", "baseline", Sys.Date())
figure_path <- here("figures", "baseline", Sys.Date())
ifelse(!dir.exists(output_path), dir.create(output_path, recursive = TRUE), FALSE)
ifelse(!dir.exists(figure_path), dir.create(figure_path, recursive = TRUE), FALSE)

### Prepare model inputs
r0s <- seq(1,5, by = 0.5)

if(eq_type == "serotype"){
  n_states <- 7522
  n_strains <- 4
  denv_mod <- odin2::odin(here("code/denv-serotype-model.R"), debug = FALSE, skip_cache = TRUE) 
  
} else if(eq_type == "strain"){
  n_states <- 27024
  n_strains <- 20
  denv_mod <- odin2::odin(here("code/denv-strain-model.R"), debug = FALSE, skip_cache = TRUE) 
  
}

### Run equilibration

for(r0 in r0s){ 
cat(paste0("Running model equilibration for r0: ", r0, "\n"))
baseline_states <- matrix(NA, nrow = n_states, ncol = 0)
for(sim in 1:50){
n_particles <- 1
calibration_years <- 100
brazil_demog <- load_demography()
demog <- wrangle_demography(brazil_demog, year_start = 1974, year_end = 2024, pad_left = max(calibration_years - 51, 0))

## Load generative model
start_time <- Sys.time()
denv_sys <- dust_system_create(denv_mod, 
                          pars = list(n_age = 81,
                                      r0 = r0,
                                      gamma = 0.2,
                                      nu = 1/365,
                                      demog = demog,
                                      scenario_years = ncol(demog),
                                      N_init = demog[,1],
                                      seasonality = 0,
                                      ext_foi = 1e-8, # external FOI to ensure strains don't go extinct during calibration
                                      amp_seas = 0.2,
                                      phase_seas = 1.56,
                                      wol_on = 0,
                                      wol_inhib = rep(1,n_strains)),
                          n_particles = n_particles)

dust_system_set_state_initial(denv_sys) # use initial conditions defined in code
t <- c(seq(1, 365*calibration_years, by = 28), 365*calibration_years) 
denv_out <- dust_system_simulate(denv_sys, t)
denv_out <- dust_unpack_state(denv_sys, denv_out)

denv_final <- dust_system_simulate(denv_sys, 365*calibration_years)
baseline_states <- cbind(baseline_states, denv_final)
qsave(baseline_states, here(output_path, paste0(eq_type, "_baseline-states_r0-", r0, ".qs")))
print(Sys.time() - start_time)
cat(paste0("Finished simulation: ", sim, "\n"))

## Wrangle outputs for plotting

inf_out <- as.data.frame(t(denv_out$inf))  |> 
  mutate(time = t)  |> 
  pivot_longer(
    cols = -time,
    names_to = "strain",
    values_to = "value")  |> 
  mutate(strain = as.numeric(gsub("^V", "", strain))) |> 
  group_by(time) |> 
  mutate(total_inf = sum(value), proportion = value/total_inf) |> 
  ungroup() 

immune_out <- as.data.frame(t(denv_out$prior_infection)) |> 
  mutate(time = t, time_step = row_number()) |> 
  pivot_longer(
    cols = -c(time, time_step),
    names_to = "prior_infection",
    values_to = "value") |> 
  mutate(prior_infection = as.numeric(gsub("^V", "", prior_infection)) - 1) |> 
  group_by(time) |> 
  mutate(total = sum(value)) |> 
  ungroup() |> 
  mutate(proportion = value/total) 
  
serotype_out <- data.frame(prior_denv1 = denv_out$prior_denv1, prior_denv2 = denv_out$prior_denv2, 
                           prior_denv3 = denv_out$prior_denv3, prior_denv4 = denv_out$prior_denv4) |>
  mutate(time = t, total_population = colSums(denv_out$N)) |>  
  pivot_longer(cols = c("prior_denv1", "prior_denv2", "prior_denv3", "prior_denv4")) |> 
  mutate(proportion = value/total_population)

## Check outputs

# Strains over time
strain_plot <- ggplot(inf_out) + 
  geom_line(aes(x = time, y = proportion, col = as.factor(strain), group = as.factor(strain))) +
  labs(x = "Time", y = "Strain prevalence") +
  theme(legend.position = "none")

infection_plot <- ggplot(inf_out) + 
  geom_line(aes(x = time, y = value, col = as.factor(strain), group = as.factor(strain))) +
  labs(x = "Time", y = "Infections") +
  theme(legend.position = "none")

# Immunity over time
immune_plot <- ggplot(immune_out) +
  geom_area(aes(x = time, y = proportion,  fill = as.factor(prior_infection), group = as.factor(prior_infection))) +
  labs(x = "Time", y = "Proportion", fill = "Number of prior infections")

# Serotype over time
serotype_plot <- ggplot(serotype_out) +
  geom_line(aes(x = time, y = proportion, col = name, group = name)) +
  labs(x = "Time", y = "Proportion of the population", color = "Prior serotype immunity")


baseline_plots <- ggarrange(infection_plot, strain_plot, immune_plot, serotype_plot, nrow = 4)
baseline_plots <- annotate_figure(baseline_plots, top= text_grob(paste0("Baseline plots with R0 = ", r0), family = plot_font))
ggsave(here(figure_path, paste0(eq_type,"_baseline-plots_r0-", r0, "-sim_", sim, ".png")), width = 180, height = 240, unit = "mm", dpi = 300, plot = baseline_plots)
rm(denv_out, denv_final, denv_sys)
}
}
