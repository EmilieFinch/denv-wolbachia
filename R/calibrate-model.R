# Script to run denv simulations

source(here("R", "utils.R"))
set.seed = 84
# Prepare model inputs

brazil_demog <- load_demography()
demog <- wrangle_demography(brazil_demog, year_start = 1974, year_end = 2024, pad_left = 49)
n_particles <- 1

## Load generative model
denv_mod <- odin2::odin(here("R/denv-strain-model.R"), debug = TRUE, skip_cache = TRUE) 
start_time <- Sys.time()
denv_sys <- dust_system_create(denv_mod, 
                          pars = list(n_age = 81,
                                      r0 = 1.5,
                                      gamma = 0.2,
                                      nu = 1/365,
                                      demog = demog,
                                      scenario_years = ncol(demog),
                                      N_init = demog[,1],
                                      ext_foi = 1e-8, # external FOI to ensure strains don't go extinct during calibration
                                      amp_seas = 0.2,
                                      phase_seas = 1.56,
                                      wol_on = 0,
                                      wol_start = 0,
                                      wol_inhib = rep(0,80)),
                          n_particles = n_particles,
                          seed = 84)

dust_system_set_state_initial(denv_sys) # use initial conditions defined in code
t <- seq(1, 365*100, by = 28) #365*100 
denv_out <- dust_system_simulate(denv_sys, t)
denv_out <- dust_unpack_state(denv_sys, denv_out)

print(Sys.time() - start_time)

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
  left_join(data.frame(total_population = colSums(denv_out$N), time_step = 1:dim(denv_out$N)[2])) |>
  mutate(proportion = value/total_population) |> 
  ungroup() 

serotype_out <- data.frame(prior_denv1 = denv_out$prior_denv1, prior_denv2 = denv_out$prior_denv2, 
                           prior_denv3 = denv_out$prior_denv3, prior_denv4 = denv_out$prior_denv4) |> 
  mutate(time = t, total_population = colSums(denv_out$N)) |> 
  pivot_longer(cols = c("prior_denv1", "prior_denv2", "prior_denv3", "prior_denv4")) |> 
  mutate(proportion = value/total_population)

## Check outputs

## Strains over time
ggplot(inf_out) + 
  geom_line(aes(x = time, y = proportion, col = as.factor(strain), group = as.factor(strain))) +
  labs(x = "Time", y = "Infections") +
  theme(legend.position = "none")

## Immunity over time

ggplot(immune_out) +
  geom_area(aes(x = time, y = proportion,  fill = as.factor(prior_infection), group = as.factor(prior_infection))) +
  labs(x = "Time", y = "Proportion", fill = "Number of prior infections")

## Serotype over time

ggplot(serotype_out) +
  geom_line(aes(x = time, y = proportion, col = name, group = name))

