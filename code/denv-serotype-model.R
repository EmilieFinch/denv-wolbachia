#### DENV serotype discrete time stochastic compartmental model ####

### Structure of script
## 1. Model parameters
## 2. Core equations
## 3. Compute probabilities of transition
## 4. Draw numbers moving between compartments
## 5. Ageing & mortality
## 6. Burden tracking
## 7. Initial states and dimensions

### Compartments are stratified by 
## i for the age group (number of age groups specified by n_age)
## Compartments C , I and R are also stratified by 
## j for immune history
## The levels of this are:
## 1, 2, 3, 4, 12, 13, 14, 23, 24, 34, 123, 124, 134, 234, 1234
## For I these are different immune history can be 0 and no one with immune history 1234 can be reinfected
## 0, 1, 2, 3, 4, 12, 13, 14, 23, 24, 34, 123, 124, 134, 234
# Note that we don't keep track of the order of past infections (so 1 then 2 then 3 is stored in the same history compartment as 3 then 2 then 1 etc)
## Compartment I is also stratified by
## k for current infecting serotypes

#### User-input model parameters ####
## Time-keeping
initial(sim_year) <- 1
update(sim_year) <- floor((time+1)/365) + 1 # track year of simulation, assuming all years are 365 days long. Note as first time is 0, sim_year == 2 happens when time = 365
scenario_years <- parameter() # number of years to simulate

## Demog & serotype parameters
n_age <- parameter(80)
n_histories <- 15
n_strains <- 4 # total number of strains modelled
r0 <- parameter() # input strain-specific R0 (e.g. each strain has this R0 when introduced alone)
gamma <- parameter(0.2)
nu <- parameter(0.002739726)
N_init <- parameter() # initial population size by age
demog <- parameter() # array of population size [age_group, sim_year]
ageing_day <- (time %% 365 == 0 && time != 0) # population ages once a year
initial(day_of_year, zero_every = 365) <- 0
update(day_of_year) <- day_of_year + 1
ext_foi <- parameter()
seasonality <- parameter(0)
amp_seas <- parameter(0.2)
phase_seas <- parameter(1.56)

## Wolbachia parameters
wol_on <- parameter() # parameter indicating whether a wolbachia intervention is being implemented
wol_inhib <- parameter() # level of wolbachia inhibition (by strain)

#### Core equations ####
update(S[]) <- floor((if(i == 1 && ageing_day) births[day_of_year] else 
  if(i == 1) births[day_of_year] + S[i] - S_out[i] else
  if(ageing_day && i > 1) S[i-1] - S_out[i-1] else
  S[i] - S_out[i]) * shift[i])

update(I[,1,]) <- floor((if(ageing_day && i == 1) 0 else
  if(ageing_day && i > 1)  I[i-1,j,k] - I_out[i-1,j,k] + n_SI[i-1,k] else 
    I[i,j,k] - I_out[i,j,k] + n_SI[i,k]) * shift[i])

update(I[,2:n_histories,]) <- floor((if(ageing_day && i == 1) 0 else
  if(ageing_day && i > 1) I[i-1,j,k] - I_out[i-1,j,k] + n_RI[i-1,j-1,k] else
    I[i,j,k] - I_out[i,j,k] + n_RI[i,j-1,k]) *shift[i])

update(C[,]) <- floor((if(ageing_day && i == 1) 0 else
 if(ageing_day && i > 1) C[i-1,j] - n_CR[i-1,j] + n_IC[i-1,j] else
    C[i,j] - n_CR[i,j] + n_IC[i,j])*shift[i])

update(R[,]) <- floor((if(ageing_day && i == 1) 0 else
  if(ageing_day && i > 1) R[i-1,j] - R_out[i-1,j] + n_CR[i-1,j] else
    R[i,j] - R_out[i,j] + n_CR[i,j])*shift[i])

update(N[]) <- if(ageing_day) demog[i, sim_year] else N[i] # don't use sum of S, I, C and R to avoid introducing lag in N

#### Calculate individual probabilities of transition ####
beta[] <- r0*gamma
beta_adj[] <- if(wol_on == 1) wol_inhib[i]*beta[i] else beta[i]
update(r0_out[]) <- beta_adj[i]/gamma

lambda_local[] <- beta_adj[i] * sum(I[,,i])/sum(N[])  # strain specific lambda

lambda[] <- if(seasonality == 1) lambda_local[i] * (1 + amp_seas * cos(2 * 3.14159 * time/365 - phase_seas)) + ext_foi else lambda_local[i] + ext_foi # add seasonality to lambda
lambda_total <- sum(lambda[])

p_IC <- 1 - exp(-gamma*dt) # I to C
p_CR <- 1 - exp(-nu*dt) # C to R

print("sim_year: {sim_year}}", when = time %% 365 == 0)

#### Draw number moving between compartments ####
## Note chain binomials for n_SI and n_RI at the end of the script 
S_out[] <- Binomial(S[i], 1 - exp(-lambda_total*dt)) 
I_out[,,] <- Binomial(I[i,j,k], p_IC)

n_IC[, 1] <- Binomial(sum(I_out[i, 1, 1]), p_IC) # seronaive infected with DENV1
n_IC[, 2] <- Binomial(sum(I_out[i, 1, 2]), p_IC) # seronaive infected with DENV2
n_IC[, 3] <- Binomial(sum(I_out[i, 1, 3]), p_IC) # seronaive infected with DENV3
n_IC[, 4] <- Binomial(sum(I_out[i, 1, 4]), p_IC) # seronaive infected with DENV4
n_IC[, 5] <- Binomial(sum(I_out[i, 2, 2]) + sum(I_out[i, 3, 1]), p_IC) # Prev DENV1 infected with DENV2 and vice versa
n_IC[, 6] <- Binomial(sum(I_out[i, 2, 3]) + sum(I_out[i, 4, 1]), p_IC) # Prev DENV1 infected with DENV3 and vice versa
n_IC[, 7] <- Binomial(sum(I_out[i, 2, 4]) + sum(I_out[i, 5, 1]), p_IC) # Prev DENV1 infected with DENV4 and vice versa
n_IC[, 8] <- Binomial(sum(I_out[i, 3, 3]) + sum(I_out[i, 4, 2]), p_IC) # Prev DENV2 infected with DENV3 and vice versa
n_IC[, 9] <- Binomial(sum(I_out[i, 3, 4]) + sum(I_out[i, 5, 2]), p_IC) # Prev DENV2 infected with DENV4 and vice versa
n_IC[, 10] <- Binomial(sum(I_out[i, 4, 4]) + sum(I_out[i, 5, 3]), p_IC) # Prev DENV3 infected with DENV4 and vice versa
n_IC[, 11] <- Binomial(sum(I_out[i, 6, 3]) + sum(I_out[i, 7, 2]) + sum(I_out[i, 9, 1]), p_IC) # Prev 2&3 infected with 1, Prev 1&2 infected with 3 and prev 1&3 infected with 2
n_IC[, 12] <- Binomial(sum(I_out[i, 6, 4]) + sum(I_out[i, 8, 2]) + sum(I_out[i, 10, 1]), p_IC) # Prev 2&4 infected with 1, Prev 1&4 infected with 2 and prev 1&2 infected with 4
n_IC[, 13] <- Binomial(sum(I_out[i, 7, 4]) + sum(I_out[i, 8, 3]) + sum(I_out[i, 11, 1]), p_IC) # Prev 3&4 infected with 1, Prev 1&4 infected with 3 and prev 1&3 infected with 4
n_IC[, 14] <- Binomial(sum(I_out[i, 9, 4]) + sum(I_out[i, 10, 3]) + sum(I_out[i, 11, 2]), p_IC) # Prev 3&4 infected with 2, Prev 2&4 infected with 3 and prev 2&3 infected with 4
n_IC[, 15] <- Binomial(sum(I_out[i, 12, 4]) + sum(I_out[i, 13, 3]) + sum(I_out[i, 14, 2]) + sum(I_out[i, 15, 1]), p_IC) # Prev 1&2&3 infected with 4, prev 1&3&4 infected with 2, prev 2&3&4 infected with 1
n_CR[,] <- Binomial(C[i,j], p_CR)

R_out[, 1] <- Binomial(R[i, 1], 1 - exp(- (lambda[2] + lambda[3] + lambda[4])*dt))
R_out[, 2] <- Binomial(R[i, 2], 1 - exp(- (lambda[1] + lambda[3] + lambda[4])*dt))
R_out[, 3] <- Binomial(R[i, 3], 1 - exp(- (lambda[1] + lambda[2] + lambda[4])*dt))
R_out[, 4] <- Binomial(R[i, 4], 1 - exp(- (lambda[1] + lambda[2] + lambda[3])*dt))
R_out[, 5] <- Binomial(R[i, 5], 1 - exp(- (lambda[3] + lambda[4])*dt))
R_out[, 6] <- Binomial(R[i, 6], 1 - exp(- (lambda[2] + lambda[4])*dt))
R_out[, 7] <- Binomial(R[i, 7], 1 - exp(- (lambda[2] + lambda[3])*dt))
R_out[, 8] <- Binomial(R[i, 8], 1 - exp(- (lambda[1] + lambda[4])*dt))
R_out[, 9] <- Binomial(R[i, 9], 1 - exp(- (lambda[1] + lambda[3])*dt))
R_out[, 10] <- Binomial(R[i, 10], 1 - exp(- (lambda[1] + lambda[2])*dt))
R_out[, 11] <- Binomial(R[i, 11], 1 - exp(- (lambda[4])*dt))
R_out[, 12] <- Binomial(R[i, 12], 1 - exp(- (lambda[3])*dt))
R_out[, 13] <- Binomial(R[i, 13], 1 - exp(- (lambda[2])*dt))
R_out[, 14] <- Binomial(R[i, 14], 1 - exp(- (lambda[1])*dt)) # chain binomial for n_RI below

#### Ageing & demography #### 

births_per_day <- demog[1,sim_year] / 365 ## spread births out throughout the year
remainder <- demog[1,sim_year] %% 365
births[] <- if(sim_year != 1 && i <= remainder) births_per_day + 1 else if(sim_year!= 1 && i > remainder) births_per_day else 0

## each year adjust each compartment to match demography data
current_N[] <- S[i] + sum(I[i,,]) + sum(C[i,]) + sum(R[i,])
shift[2:n_age] <- if(ageing_day && current_N[i-1] > 0) demog[i, sim_year]/current_N[i-1] else 1.0 # adjust population sizes so they match demographic data from upcoming year
shift[1] <- 1 # for youngest age group all enter as susceptible

#### Track burden outputs ####
## Infections by strain

update(inf[]) <- sum(I[,,i])
update(inf_age[]) <- sum(I[i,,])

## Immune states
update(prior_infection[1]) <- sum(S[]) # sum seronaitve
update(prior_infection[2]) <- sum(R[,1:4]) # sum one previous infection
update(prior_infection[3]) <- sum(R[,5:10]) # sum two previous infections
update(prior_infection[4]) <- sum(R[,11:14]) # sum three previous infections
update(prior_infection[5]) <- sum(R[,15]) # sum four previous infections

update(prior_denv1) <-sum(R[, 1]) + sum(R[, 5]) + sum(R[, 6]) + sum(R[, 7]) + sum(R[, 11]) + sum(R[, 12]) + sum(R[, 13]) + sum(R[, 15]) # sum all those with immunity to DENV-1
update(prior_denv2) <- sum(R[, 2]) + sum(R[, 5]) + sum(R[, 8]) + sum(R[, 9]) + sum(R[, 11]) + sum(R[, 12]) + sum(R[, 14]) + sum(R[, 15]) # sum all those with immunity to DENV-2
update(prior_denv3) <- sum(R[, 3]) + sum(R[, 6]) + sum(R[, 8]) + sum(R[, 10]) + sum(R[, 11]) + sum(R[, 13]) + sum(R[, 14]) + sum(R[, 15]) # sum all those with immunity to DENV-3
update(prior_denv4) <- sum(R[, 4]) + sum(R[, 7]) + sum(R[, 9]) + sum(R[, 10]) + sum(R[, 12]) + sum(R[, 13]) + sum(R[, 14]) + sum(R[, 15]) # sum all those with immunity to DENV-4


#### Initial states & dimensions ####
initial(N[]) <- N_init[i]
initial(S[]) <- if(i == 11) N_init[i] - 4 else N_init[i]
initial(I[11,1,1]) <- 1
initial(I[11,2,2]) <- 1
initial(I[11,3,3]) <- 1
initial(I[11,4,4]) <- 1
initial(C[,]) <- 0
initial(R[,]) <- 0
initial(r0_out[]) <- r0
initial(inf[]) <- 0
initial(inf_age[]) <- 0
initial(prior_infection[]) <- 0
initial(prior_denv1) <- 0
initial(prior_denv2) <- 0
initial(prior_denv3) <- 0
initial(prior_denv4) <- 0

dim(S) <- n_age
dim(I) <- c(n_age, n_histories, n_strains)
dim(C) <- c(n_age, n_histories)
dim(R) <- c(n_age, n_histories)
dim(S_out) <- n_age
dim(N_init) <- n_age
dim(beta) <- n_strains
dim(beta_adj) <- n_strains
dim(r0_out) <- n_strains
dim(wol_inhib) <- n_strains
dim(lambda_local) <- n_strains
dim(lambda) <- n_strains
dim(I_out) <- c(n_age, n_histories, n_strains)
dim(R_out) <- c(n_age, n_histories)
dim(n_SI) <- c(n_age, n_strains)
dim(n_RI) <- c(n_age, n_histories, n_strains)
dim(n_IC) <- c(n_age, n_histories)
dim(n_CR) <- c(n_age, n_histories)
dim(N) <- n_age
dim(inf) <- n_strains
dim(inf_age) <- n_age
dim(prior_infection) <- c(5)
dim(demog) <- c(n_age, scenario_years)
dim(births) <- c(365)
dim(shift) <- c(n_age)
dim(current_N) <- c(n_age)

## Chain binomial for strain specific infections ## 
## From susceptible
n_SI[, 1] <- if (sum(lambda[1:4]) > 0) Binomial(S_out[i], lambda[1] / sum(lambda[1:4])) else 0
n_SI[, 2] <- if (sum(lambda[2:4]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1]), lambda[2] / sum(lambda[2:4])) else 0
n_SI[, 3] <- if (sum(lambda[3:4]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:2]), lambda[3] / sum(lambda[3:4])) else 0
n_SI[, 4] <- max(S_out[i] - sum(n_SI[i, 1:3]),0)


n_RI[, 1, 2] <- if ((sum(lambda[2:4])) > 0) Binomial(R_out[i, 1], lambda[2] / (sum(lambda[2:4]))) else 0
n_RI[, 1, 3] <- if ((sum(lambda[3:4])) > 0) Binomial(R_out[i, 1] - (n_RI[i, 1, 2]), lambda[3] / (sum(lambda[3:4]))) else 0
n_RI[, 1, 4] <- R_out[i, 1] - (sum(n_RI[i, 1, 2:3]))
n_RI[, 2, 1] <- if ((lambda[1] + sum(lambda[3:4])) > 0) Binomial(R_out[i, 2], lambda[1] / (lambda[1] + sum(lambda[3:4]))) else 0
n_RI[, 2, 3] <- if ((sum(lambda[3:4])) > 0) Binomial(R_out[i, 2] - (n_RI[i, 2, 1]), lambda[3] / (sum(lambda[3:4]))) else 0
n_RI[, 2, 4] <- R_out[i, 2] - (n_RI[i, 2, 1] + n_RI[i, 2, 3])
n_RI[, 3, 1] <- if ((sum(lambda[1:2]) + lambda[4]) > 0) Binomial(R_out[i, 3], lambda[1] / (sum(lambda[1:2]) + lambda[4])) else 0
n_RI[, 3, 2] <- if ((lambda[2] + lambda[4]) > 0) Binomial(R_out[i, 3] - (n_RI[i, 3, 1]), lambda[2] / (lambda[2] + lambda[4])) else 0
n_RI[, 3, 4] <- R_out[i, 3] - (sum(n_RI[i, 3, 1:2]))
n_RI[, 4, 1] <- if ((sum(lambda[1:3])) > 0) Binomial(R_out[i, 4], lambda[1] / (sum(lambda[1:3]))) else 0
n_RI[, 4, 2] <- if ((sum(lambda[2:3])) > 0) Binomial(R_out[i, 4] - (n_RI[i, 4, 1]), lambda[2] / (sum(lambda[2:3]))) else 0
n_RI[, 4, 3] <- R_out[i, 4] - (sum(n_RI[i, 4, 1:2]))
n_RI[, 5, 3] <- if ((sum(lambda[3:4])) > 0) Binomial(R_out[i, 5], lambda[3] / (sum(lambda[3:4]))) else 0
n_RI[, 5, 4] <- R_out[i, 5] - (n_RI[i, 5, 3])
n_RI[, 6, 2] <- if ((lambda[2] + lambda[4]) > 0) Binomial(R_out[i, 6], lambda[2] / (lambda[2] + lambda[4])) else 0
n_RI[, 6, 4] <- R_out[i, 6] - (n_RI[i, 6, 2])
n_RI[, 7, 2] <- if ((sum(lambda[2:3])) > 0) Binomial(R_out[i, 7], lambda[2] / (sum(lambda[2:3]))) else 0
n_RI[, 7, 3] <- R_out[i, 7] - (n_RI[i, 7, 2])
n_RI[, 8, 1] <- if ((lambda[1] + lambda[4]) > 0) Binomial(R_out[i, 8], lambda[1] / (lambda[1] + lambda[4])) else 0
n_RI[, 8, 4] <- R_out[i, 8] - (n_RI[i, 8, 1])
n_RI[, 9, 1] <- if ((lambda[1] + lambda[3]) > 0) Binomial(R_out[i, 9], lambda[1] / (lambda[1] + lambda[3])) else 0
n_RI[, 9, 3] <- R_out[i, 9] - (n_RI[i, 9, 1])
n_RI[, 10, 1] <- if ((sum(lambda[1:2])) > 0) Binomial(R_out[i, 10], lambda[1] / (sum(lambda[1:2]))) else 0
n_RI[, 10, 2] <- R_out[i, 10] - (n_RI[i, 10, 1])
n_RI[, 11, 4] <- R_out[i,11] 
n_RI[, 12, 3] <- R_out[i,12] 
n_RI[, 13, 2] <- R_out[i,13] 
n_RI[, 14, 1] <- R_out[i, 14]