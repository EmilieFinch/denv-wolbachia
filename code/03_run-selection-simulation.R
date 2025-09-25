#### Code to run model simulations for strain selection analysis ####

### This script runs simulations over 75 years assuming different baseline R0 values using denv-strain-model in odin
## Inputs:
# The model uses equilibrated states obtained from 01_equilibrate-model. R as initial conditions (with a different set of initial conditions for each stochastic simulation)
# User inputs: r0 level and Wolbachia strain (wmel or walbb)
## We then draw 5 strains are drawn for each serotype from the dataset along with their corresponding relative dissemination under the Wolbachia strain (resulting in 20 strains in total within the model)
## We then simulate strain dynamics over 75 years following Wolbachia introduction 
## Outputs saved are: 1) quantiles of simulation outputs, 2) a csv of the proportion of infections belonging to each strain over the simulation period

# Note that this script is set up to run on an HPC system as an array job, with the inhibition level passed as a command line argument using the csv files
# strain-specific-input-values.csv and reemergence-input-values.csv

library(here)
source(here("code", "00_load-packages.R"))
source(here("code", "utils.R"))
output_path <- here("output", "simulation", Sys.Date())
ifelse(!dir.exists(output_path), dir.create(output_path, recursive = TRUE), FALSE)

#### config
set.seed(27)
#c_args = commandArgs(trailingOnly = TRUE)

r0 <- 2 #as.numeric(c_args[[1]])
mosquito <- "wmel" #as.character(c_args[[2]])
n_particles <- 50
simulation_years <- 75

#### Load model
denv_mod <- odin2::odin(here("code/denv-strain-model.R"), debug = TRUE, skip_cache = TRUE)

#### Prepare model inputs
brazil_demog <- load_demography()
demog <- wrangle_demography(brazil_demog, year_start = 2025, year_end = 2100, pad_left = 0)
initial_states <- qread(here("output", "baseline", "2025-09-10", paste0("strain_baseline-states_r0-", r0, ".qs")))

inhib_table <- qread(here("data", paste0(mosquito, "-inhib.qs")))

for(draw in 1:20){
  
  inhib_selected <- inhib_table |>
    group_by(serotype) |>
    sample_n(size = 5, replace = FALSE)
  
  inhib_values <- inhib_selected |>
    pull(dissemination_norm)
  
  ### Set up model
  start_time <- Sys.time()
  denv_sys <- dust_system_create(denv_mod,
                                 pars = list(n_age = 81,
                                             r0 = r0,
                                             gamma = 0.2,
                                             nu = 1/365,
                                             demog = demog,
                                             scenario_years = ncol(demog),
                                             N_init = demog[,1],
                                             seasonality = 1,
                                             ext_foi = 1e-8, # external FOI set to 0 for simulation
                                             amp_seas = 0.2,
                                             phase_seas = 1.56,
                                             wol_on = 1,
                                             wol_inhib = inhib_values),
                                 seed = 27,
                                 n_particles = n_particles,
                                 n_threads = 10)
  
  idx <- c("prior_infection" = dust_unpack_index(denv_sys)$prior_infection, "prior_denv1" = dust_unpack_index(denv_sys)$prior_denv1,
           "prior_denv2" = dust_unpack_index(denv_sys)$prior_denv2, "prior_denv3" = dust_unpack_index(denv_sys)$prior_denv3, "prior_denv4" = dust_unpack_index(denv_sys)$prior_denv4,
           "inf" = dust_unpack_index(denv_sys)$inf) # set indices to return
  dust_system_set_state(denv_sys, initial_states) # use initial conditions defined in code
  t <- seq(1, 365*simulation_years, by = 1)
  
  cat(paste0("Running strain simulations"))
  sim_raw <- dust_system_simulate(denv_sys, t, index_state = idx)
  
  ## Wrangle and save output
  cat(paste0("Saving simulation output"))
  
  # Raw simulations
  sim_out <- as.data.frame.table(sim_raw) |>
    rename(state = Var1, run = Var2, time = Var3, value = Freq) |>
    mutate(run = as.numeric(run), time = as.numeric(time)) |>
    extract(state, into = c("name", "level"), regex = "([a-zA-Z]+)([0-9]+)", remove = FALSE) |>
    mutate(draw = draw, mosquito = mosquito)
  
  # Quantiles
  quantiles_out <- sim_out |>
    group_by(state, name, level, time, draw, mosquito) |>
    summarise(median = median(value), q_0.25 = quantile(value, 0.25), q_0.75 = quantile(value, 0.75),
              q_0.025 = quantile(value, 0.025), q_0.975 = quantile(value, 0.975)) |> 
    mutate(year =  floor((time+1)/365) + 2025) |>
    left_join(brazil_demog |> group_by(year) |> summarise(population = sum(population))) 
    
  qsave(quantiles_out, here(output_path, paste0("strain-quantile-out_", mosquito, "_", draw,".qs")))
  
  plot_times <-c(1, as.vector(outer(c(182, 365), 0:74 * 365, "+")))
 
  # Proportions of infections by strain
  inf_proportions <- sim_out |>
    filter(name == "inf") |>
    filter(time %in% plot_times) |>
    group_by(name, run, time) |>
    mutate(total_infections = sum(value)) |>
    group_by(name, run, time, level, total_infections) |>
    mutate(proportions = value/total_infections) |>
    mutate(virus = unique(inhib_selected$virus)[as.numeric(level)]) |>
    mutate(serotype = case_when(as.numeric(level) <= 5 ~ 1,
                                as.numeric(level) <= 10 ~ 2,
                                as.numeric(level) <= 15 ~ 3,
                                T ~ 4)) |>
    mutate(draw = draw, mosquito = mosquito)
  
  
  if (file.exists(here(output_path, "strain-out_", mosquito, ".csv"))) {
    write_csv(inf_proportions, here(output_path, paste0("strain-out_", mosquito, ".csv")),append = TRUE)
  } else {
    write_csv(inf_proportions, here(output_path, paste0("strain-out_", mosquito, ".csv")))
  }
  
  
}

  