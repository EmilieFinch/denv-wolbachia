## Run model simulations
library(here)
source(here("code", "00_load-packages.R"))
source(here("code", "utils.R"))
output_path <- here("output", "simulation", Sys.Date())
ifelse(!dir.exists(output_path), dir.create(output_path, recursive = TRUE), FALSE)
figure_path <- here("figures", "simulation", Sys.Date())
ifelse(!dir.exists(figure_path), dir.create(figure_path, recursive = TRUE), FALSE)

# config
n_particles <- 50
simulation_years <- 75
c_args = commandArgs(trailingOnly = TRUE)

r0 <-as.numeric(c_args[[1]])
inhibs <- seq(0.9,0, by = -0.1)

# load model
denv_mod <- odin2::odin(here("code/denv-serotype-model.R"), debug = TRUE, skip_cache = TRUE)

for(inhib_level in inhib_levels){ 
for(serotype_varying in 1:4){
  # Prepare model inputs
  brazil_demog <- load_demography()
  demog <- wrangle_demography(brazil_demog, year_start = 2025, year_end = 2100, pad_left = 0)
  
  initial_states <- qread(here("output", "calibration", "2025-06-06", paste0("calibration-states_r0-", r0, ".qs")))
  
  wol_inhib <- rep(0,4)
  wol_inhib[serotype_varying] <- inhib_level
  
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
                                             wol_inhib = wol_inhib),
                                 n_particles = n_particles,
                                 n_threads = 10)
  
  idx <- c("prior_infection" = dust_unpack_index(denv_sys)$prior_infection, "prior_denv1" = dust_unpack_index(denv_sys)$prior_denv1,
           "prior_denv2" = dust_unpack_index(denv_sys)$prior_denv2, "prior_denv3" = dust_unpack_index(denv_sys)$prior_denv3, "prior_denv4" = dust_unpack_index(denv_sys)$prior_denv4,
           "inf" = dust_unpack_index(denv_sys)$inf) # set indices to return
  dust_system_set_state(denv_sys, as.vector(initial_states)) # use initial conditions defined in code
  t <- seq(1, 365*simulation_years, by = 1)
  
  cat(paste0("Running simulations for r0: ", r0, " and inhibition level: ", inhib_level, " and serotype: ", serotype_varying, "\n"))
  sim_raw <- dust_system_simulate(denv_sys, t, index_state = idx)
  
  ## Wrangle and save output
  cat(paste0("Saving output for r0: ", r0, " and inhibition level: ", inhib_level, " and serotype: ", serotype_varying, "\n"))
  
  ### Raw simulations
  sim_out <- as.data.frame.table(sim_raw) |>
    rename(state = Var1, run = Var2, time = Var3, value = Freq) |>
    mutate(run = as.numeric(run), time = as.numeric(time)) |>
    extract(state, into = c("name", "level"), regex = "([a-zA-Z]+)([0-9]+)", remove = FALSE) |>
    mutate(inhib_level = inhib_level, r0 = r0, serotype_varying = serotype_varying)
  
  qsave(sim_out,here(output_path, paste0("simulation-out_r0-",r0, "_inhib-level-", inhib_level,"_serotype-varying-", serotype_varying, ".qs")))
  
  ### Time to emergence
  time_to_emergence <- sim_out |>
    filter(name == "inf") |>
    group_by(name, run, time, inhib_level, r0, serotype_varying) |>
    filter(level == serotype_varying) |>
    mutate(year =  floor((time+1)/365) + 2025) |>
    left_join(brazil_demog |> group_by(year) |> summarise(population = sum(population))) |>
    mutate(inc = (value/population)*100000) |>
    group_by(run) |>
    mutate(flag = inc > 1,
           extinction = case_when(flag == FALSE & flag != lag(flag)  ~ T,
                                  T ~ F),
           extinction = if_else(row_number() == which(extinction)[1], 1L, 0L),
           reemergence = case_when(flag == TRUE & flag != lag(flag) ~ T,
                                   T ~ F),
           reemergence = if_else(row_number() == which(reemergence)[1], 1L, 0L),
           time_to_extinction = time[extinction == 1],
           time_to_reemergence = time[reemergence == 1])
  
  emergence_out <- time_to_emergence |>
    mutate(date = as.Date("2025-01-01") + time - 1) |>
    select(name, level, run, inhib_level, r0,serotype_varying, time_to_extinction, time_to_reemergence) |>
    unique() |>
    arrange(run)
  
  if (file.exists(here(output_path, paste0("time-to-emergence_r0-",r0, ".csv")))) {
    write_csv(emergence_out, here(output_path, paste0("time-to-emergence_r0-",r0, ".csv")),append = TRUE)
  } else {
    write_csv(emergence_out, here(output_path, paste0("time-to-emergence_r0-",r0, ".csv")))
  }
  
  
  ### Infection plot
  infection_plot <- sim_out |>
    filter(name == "inf") |>
    mutate(date = as.Date("2025-01-01") + time - 1) |>
    group_by(name, level, time, date, inhib_level, r0) |>
    summarise(median = median(value), q_0.025 = quantile(value, 0.025), q_0.975 = quantile(value, 0.975),
              q_0.25 = quantile(value, 0.25), q_0.75 = quantile(value, 0.75)) |>
    ungroup() |>
    ggplot() +
    geom_line(aes(x = date, y = median, col = as.factor(level), group = as.factor(level))) +
    geom_ribbon(aes(x = date, ymin = q_0.25, ymax = q_0.75, fill = as.factor(level), group = as.factor(level)), alpha = 0.4) +
    geom_ribbon(aes(x = date, ymin = q_0.025, ymax = q_0.975, fill = as.factor(level), group = as.factor(level)), alpha = 0.2) +
    labs(y = "Infections", x = "Date", col = "Serotype", fill = "Serotype") +
    scale_color_manual(values = c("#E07529FF", "#FAAE32FF", "#7F7991FF", "#5D4F36FF", "#B39085FF")) +
    scale_fill_manual(values = c("#E07529FF", "#FAAE32FF", "#7F7991FF", "#5D4F36FF", "#B39085FF"))
  
  ggsave(here(figure_path, paste0("infection-plot_r0-", r0, "_inhib-level-", inhib_level, "_serotype-varying-", serotype_varying, ".png")),
         width = 180, height = 150, unit = "mm", dpi = 300, plot = infection_plot)
  
  cat(Sys.time() - start_time)
  cat(paste0("Completed simulations for r0: ", r0, " and inhib level: ", inhib_level, " with serotype varying: ", serotype_varying, "\n"))
  rm(sim_raw, sim_out, quantiles_out, time_to_emergence, emergence_out, infection_plot, denv_sys)
}
}
