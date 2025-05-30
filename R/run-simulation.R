## Run model simulations

source(here("R", "utils.R"))
output_path <- here("output", "simulation")
ifelse(!dir.exists(output_path), dir.create(output_path, recursive = TRUE), FALSE)

# config
n_particles <- 50
simulation_years <- 10

r0s <- c(2, 4)
mosquitos <- c("wmel", "walbb", "wt")
for(r0 in r0s){
  for(mosquito in mosquitos){

# Prepare model inputs
brazil_demog <- load_demography()
demog <- wrangle_demography(brazil_demog, year_start = 2025, year_end = 2070, pad_left = 0)

initial_states <- qread(here("output", "calibration", paste0("calibration-states_r0-", r0, ".qs")))

if(mosquito == "wt"){
  wol_inhib <- rep(1,80)
} else {
wol_in <- qread(here("data", paste0(mosquito, "-inhib-values.qs")))
wol_inhib <- wol_in[1,]}
  
## Load generative model
denv_mod <- odin2::odin(here("R/denv-strain-model.R"), debug = TRUE, skip_cache = TRUE) 
start_time <- Sys.time()
denv_sys <- dust_system_create(denv_mod, 
                               pars = list(n_age = 81,
                                           r0 = r0,
                                           gamma = 0.2,
                                           nu = 1/365,
                                           demog = demog,
                                           scenario_years = ncol(demog),
                                           N_init = demog[,1],
                                           ext_foi = 0, # external FOI set to 0 for simulation
                                           amp_seas = 0.2,
                                           phase_seas = 1.56,
                                           wol_on = 1,
                                           wol_inhib = wol_inhib),
                               n_particles = n_particles,
                               n_threads = 4)

idx <- c("prior_infection" = dust_unpack_index(denv_sys)$prior_infection, "prior_denv1" = dust_unpack_index(denv_sys)$prior_denv1, 
         "prior_denv2" = dust_unpack_index(denv_sys)$prior_denv2, "prior_denv3" = dust_unpack_index(denv_sys)$prior_denv3, "prior_denv4" = dust_unpack_index(denv_sys)$prior_denv4, 
         "inf" = dust_unpack_index(denv_sys)$inf) # set indices to return
dust_system_set_state(denv_sys, as.vector(initial_states)) # use initial conditions defined in code
t <- seq(1, 365*simulation_years, by = 1) 

print(paste0("Running simulations for r0: ", r0, " and mosquito: ", mosquito))
sim_raw <- dust_system_simulate(denv_sys, t, index_state = idx)

## Wrangle and save output
sim_out <- as.data.frame.table(sim_raw) |> 
  rename(state = Var1, run = Var2, time = Var3, value = Freq) |> 
  mutate(run = as.numeric(run), time = as.numeric(time)) |> 
  extract(state, into = c("name", "level"), regex = "([a-zA-Z]+)([0-9]+)", remove = FALSE) 

qsave(sim_out,here(output_path, paste0("simulation-out_r0-",r0, "_", mosquito, ".qs")))

print(Sys.time() - start_time)
print(paste0("Completed simulations for r0: ", r0, " and mosquito: ", mosquito))
  }
}

