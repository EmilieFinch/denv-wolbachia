# Script to run denv simulations

# Prepare model inputs

## Load generative model
denv_mod <- odin2::odin(here("R/denv-strain-model.R"), debug = TRUE, skip_cache = TRUE) 

denv_sys <- dust_system_create(denv_mod, 
                          pars = list(n_age = 80,
                                      r0 = 2,
                                      gamma = 0.2,
                                      N_init = rep(2000000,80),
                                      I_init = I_init,
                                      wol_on = 0,
                                      wol_start = 0,
                                      wol_inhib = rep(0,60)))
dust_system_set_state_initial(denv_sys) # use initial conditions defined in code
t <- seq(1, 365*2, by = 1)
denv_out <- dust_system_simulate(denv_sys, t)
denv_out <- dust_unpack_state(denv_sys, denv_out)

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
  mutate(time = t) |> 
  pivot_longer(
    cols = -time,
    names_to = "prior_infection",
    values_to = "value") |> 
  mutate(prior_infection = as.numeric(gsub("^V", "", prior_infection)) - 1) |> 
  left_join(data.frame(total_population = colSums(denv_out$N), time = 1:dim(denv_out$N)[2])) |>
  mutate(proportion = value/total_population) |> 
  ungroup() 
  

## Check outputs

## Strains over time
ggplot(inf_out) + 
  geom_line(aes(x = time, y = value, col = as.factor(strain), group = as.factor(strain)))

## Immunity over time

ggplot(immune_out) +
  geom_area(aes(x = time, y = proportion,  fill = prior_infection, group = prior_infection))
