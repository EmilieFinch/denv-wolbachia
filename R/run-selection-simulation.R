# Run selection simulation
library(here)
source(here("code", "00_load-packages.R"))
source(here("code", "utils.R"))
output_path <- here("output", "strain-simulation", Sys.Date())
ifelse(!dir.exists(output_path), dir.create(output_path, recursive = TRUE), FALSE)
figure_path <- here("figures", "strain-simulation", Sys.Date())
ifelse(!dir.exists(figure_path), dir.create(figure_path, recursive = TRUE), FALSE)

# config
c_args = commandArgs(trailingOnly = TRUE)

r0 <-as.numeric(c_args[[1]])
mosquito <-as.character(c_args[[2]])
n_particles <- 50
simulation_years <- 75

# load model
denv_mod <- odin2::odin(here("code/denv-strain-model.R"), debug = TRUE, skip_cache = TRUE)

# Prepare model inputs
brazil_demog <- load_demography()
demog <- wrangle_demography(brazil_demog, year_start = 2025, year_end = 2100, pad_left = 0)
initial_states <- qread(here("output", "strain-calibration", "2025-06-10", paste0("calibration-states_r0-", r0, ".qs")))

inhib_table <- qread(here("data", paste0(mosquito, "-inhib.qs")))

for(draw in 1:20){
  inhib_selected <- inhib_table |>
    group_by(serotype) |>
    sample_n(size = 5, replace = FALSE)
  
  inhib_values <- inhib_selected |>
    pull(dissemination_norm)
  
  ## Set up model
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
                                 n_particles = n_particles,
                                 n_threads = 10)
  
  idx <- c("prior_infection" = dust_unpack_index(denv_sys)$prior_infection, "prior_denv1" = dust_unpack_index(denv_sys)$prior_denv1,
           "prior_denv2" = dust_unpack_index(denv_sys)$prior_denv2, "prior_denv3" = dust_unpack_index(denv_sys)$prior_denv3, "prior_denv4" = dust_unpack_index(denv_sys)$prior_denv4,
           "inf" = dust_unpack_index(denv_sys)$inf) # set indices to return
  dust_system_set_state(denv_sys, as.vector(initial_states)) # use initial conditions defined in code
  t <- seq(1, 365*simulation_years, by = 1)
  
  cat(paste0("Running strain simulations"))
  sim_raw <- dust_system_simulate(denv_sys, t, index_state = idx)
  
  ## Wrangle and save output
  cat(paste0("Saving simulation output"))
  
  ### Raw simulations
  sim_out <- as.data.frame.table(sim_raw) |>
    rename(state = Var1, run = Var2, time = Var3, value = Freq) |>
    mutate(run = as.numeric(run), time = as.numeric(time)) |>
    extract(state, into = c("name", "level"), regex = "([a-zA-Z]+)([0-9]+)", remove = FALSE) |>
    mutate(draw = draw)
  
  #qsave(sim_out,here(output_path, paste0("strain-simulation-out.qs")))
  
  ### Quantiles
  quantiles_out <- sim_out |>
    group_by(state, name, level, time,draw) |>
    summarise(median = median(value), q_0.25 = quantile(value, 0.25), q_0.75 = quantile(value, 0.75),
              q_0.025 = quantile(value, 0.025), q_0.975 = quantile(value, 0.975))
  qsave(quantiles_out, here(output_path, paste0("strain-quantile-out_", draw,".qs")))
  
  # Aim: data frame with the proportion of each of the strains in each run
  # Then plot these with colours for serotypes and total n for time periods: pre-release (calibration state), after 5 years, 20 years, 50 years 75 years
  
  plot_times <- c(1, 365*5, 365*20, 365*50, 365*75) # in days
  
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
    mutate(draw = draw)
  
  
  if (file.exists(here(output_path, "strain-proportions.csv"))) {
    write_csv(inf_proportions, here(output_path, paste0("strain-proportions.csv")),append = TRUE)
  } else {
    write_csv(inf_proportions, here(output_path, paste0("strain-proportions.csv")))
  }
  
  
}

  