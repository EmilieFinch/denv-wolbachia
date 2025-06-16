#### Definition of discrete time stochastic compartmental model
### Structure of script
## 1. Model parameters
## 2. Core equations
## 3. Compute probabilities of transition
## 4. Draw numbers moving between compartments
## 5. Ageing & mortality
## 6. Burden tracking
## 7. Initial states and dimensions

### Compartments are stratified by 
## i for the age group (number of age groups specified by n_age, with age group width specified in size_age)
## Compartments C , I and R are also stratified by 
## j for immune history
## Compartment I is also stratified by
## k for immune history
## The levels of this are:
## 1, 2, 3, 4, 12, 13, 14, 23, 24, 34, 123, 124, 134, 234, 1234
## For I these are different immune history can be 0 and no one with immune history 1234 can be reinfected
## 0, 1, 2, 3, 4, 12, 13, 14, 23, 24, 34, 123, 124, 134, 234
# Note that we don't keep track of the order of past infections


#### User-input model parameters ####
## Time-keeping
initial(sim_year) <- 1
update(sim_year) <- floor((time+1)/365) + 1 # track year of simulation, assuming all years are 365 days long. Note as first time is 0, sim_year == 2 happens when time = 365
scenario_years <- parameter() # number of years to simulate

## Demog & serotype parameters
n_age <- parameter(80)
n_histories <- 15
n_strains <- 20 #total number of strains modelled
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
wol_on <- parameter()
wol_inhib <- parameter() # level of wolbachia inhibition (by strain)

#### Core equations ####
update(S[]) <- floor((if(i == 1 && ageing_day) births[day_of_year + 1] else 
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
lambda_1 <- sum(lambda[1:5]) # DENV1 lambda
lambda_2 <- sum(lambda[6:10]) # DENV2 lambda
lambda_3 <- sum(lambda[11:15]) # DENV3 lambda
lambda_4 <- sum(lambda[16:20]) # DENV4 lambda

lambda_total <- sum(lambda[])

p_IC <- 1 - exp(-gamma*dt) # I to C
p_CR <- 1 - exp(-nu*dt) # C to R

print("sim_year: {sim_year}, r0_out[12]: {r0_out[12]}, lambda[12]: {lambda[12]}, lambda_total: {lambda_total}, I_example: {I[12,7,20]}")

#### Draw number moving between compartments ####
## Note chain binomials for n_SI and n_RI at the end of the script 
S_out[] <- Binomial(S[i], 1 - exp(-lambda_total*dt)) 
I_out[,,] <- Binomial(I[i,j,k], p_IC)

n_IC[, 1] <- Binomial(sum(I_out[i, 1, 1:5]), p_IC)
n_IC[, 2] <- Binomial(sum(I_out[i, 1, 6:10]), p_IC)
n_IC[, 3] <- Binomial(sum(I_out[i, 1, 11:15]), p_IC)
n_IC[, 4] <- Binomial(sum(I_out[i, 1, 16:20]), p_IC)
n_IC[, 5] <- Binomial(sum(I_out[i, 2, 6:10]) + sum(I_out[i, 3, 1:5]), p_IC)
n_IC[, 6] <- Binomial(sum(I_out[i, 2, 11:15]) + sum(I_out[i, 4, 1:5]), p_IC)
n_IC[, 7] <- Binomial(sum(I_out[i, 2, 16:20]) + sum(I_out[i, 5, 1:5]), p_IC)
n_IC[, 8] <- Binomial(sum(I_out[i, 3, 11:15]) + sum(I_out[i, 4, 6:10]), p_IC)
n_IC[, 9] <- Binomial(sum(I_out[i, 3, 16:20]) + sum(I_out[i, 5, 6:10]), p_IC)
n_IC[, 10] <- Binomial(sum(I_out[i, 4, 16:20]) + sum(I_out[i, 5, 11:15]), p_IC)
n_IC[, 11] <- Binomial(sum(I_out[i, 6, 11:15]) + sum(I_out[i, 7, 6:10]) + sum(I_out[i, 9, 1:5]), p_IC)
n_IC[, 12] <- Binomial(sum(I_out[i, 6, 16:20]) + sum(I_out[i, 8, 6:10]) + sum(I_out[i, 10, 1:5]), p_IC)
n_IC[, 13] <- Binomial(sum(I_out[i, 7, 16:20]) + sum(I_out[i, 8, 11:15]) + sum(I_out[i, 11, 1:5]), p_IC)
n_IC[, 14] <- Binomial(sum(I_out[i, 9, 16:20]) + sum(I_out[i, 10, 11:15]) + sum(I_out[i, 11, 6:10]), p_IC)
n_IC[, 15] <- Binomial(sum(I_out[i, 12, 16:20]) + sum(I_out[i, 13, 11:15]) + sum(I_out[i, 14, 6:10]) + sum(I_out[i, 15, 1:5]), p_IC)
n_CR[,] <- Binomial(C[i,j], p_CR)

R_out[, 1] <- Binomial(R[i, 1], 1 - exp(- (lambda_2 + lambda_3 + lambda_4)*dt))
R_out[, 2] <- Binomial(R[i, 2], 1 - exp(- (lambda_1 + lambda_3 + lambda_4)*dt))
R_out[, 3] <- Binomial(R[i, 3], 1 - exp(- (lambda_1 + lambda_2 + lambda_4)*dt))
R_out[, 4] <- Binomial(R[i, 4], 1 - exp(- (lambda_1 + lambda_2 + lambda_3)*dt))
R_out[, 5] <- Binomial(R[i, 5], 1 - exp(- (lambda_3 + lambda_4)*dt))
R_out[, 6] <- Binomial(R[i, 6], 1 - exp(- (lambda_2 + lambda_4)*dt))
R_out[, 7] <- Binomial(R[i, 7], 1 - exp(- (lambda_2 + lambda_3)*dt))
R_out[, 8] <- Binomial(R[i, 8], 1 - exp(- (lambda_1 + lambda_4)*dt))
R_out[, 9] <- Binomial(R[i, 9], 1 - exp(- (lambda_1 + lambda_3)*dt))
R_out[, 10] <- Binomial(R[i, 10], 1 - exp(- (lambda_1 + lambda_2)*dt))
R_out[, 11] <- Binomial(R[i, 11], 1 - exp(- (lambda_4)*dt))
R_out[, 12] <- Binomial(R[i, 12], 1 - exp(- (lambda_3)*dt))
R_out[, 13] <- Binomial(R[i, 13], 1 - exp(- (lambda_2)*dt))
R_out[, 14] <- Binomial(R[i, 14], 1 - exp(- (lambda_1)*dt)) # chain binomial for n_RI below

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
update(prior_infection[1]) <- sum(S[])
update(prior_infection[2]) <- sum(R[,1:4])
update(prior_infection[3]) <- sum(R[,5:10])
update(prior_infection[4]) <- sum(R[,11:14])
update(prior_infection[5]) <- sum(R[,15])

update(prior_denv1) <-sum(R[, 1]) + sum(R[, 5]) + sum(R[, 6]) + sum(R[, 7]) + sum(R[, 11]) + sum(R[, 12]) + sum(R[, 13]) + sum(R[, 15])
update(prior_denv2) <- sum(R[, 2]) + sum(R[, 5]) + sum(R[, 8]) + sum(R[, 9]) + sum(R[, 11]) + sum(R[, 12]) + sum(R[, 14]) + sum(R[, 15])
update(prior_denv3) <- sum(R[, 3]) + sum(R[, 6]) + sum(R[, 8]) + sum(R[, 10]) + sum(R[, 11]) + sum(R[, 13]) + sum(R[, 14]) + sum(R[, 15])
update(prior_denv4) <- sum(R[, 4]) + sum(R[, 7]) + sum(R[, 9]) + sum(R[, 10]) + sum(R[, 12]) + sum(R[, 13]) + sum(R[, 14]) + sum(R[, 15])


#### Initial states & dimensions ####
initial(N[]) <- N_init[i]
initial(S[]) <- if(i == 11) N_init[i] - 20 else N_init[i]
initial(I[11,1,1:5]) <- 1
initial(I[11,2,6:10]) <- 1
initial(I[11,3,11:15]) <- 1
initial(I[11,4,16:20]) <- 1
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
n_SI[, 1] <- if (sum(lambda[1:20]) > 0) Binomial(S_out[i], lambda[1] / sum(lambda[1:20])) else 0
n_SI[, 2] <- if (sum(lambda[2:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:1]), lambda[2] / sum(lambda[2:20])) else 0
n_SI[, 3] <- if (sum(lambda[3:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:2]), lambda[3] / sum(lambda[3:20])) else 0
n_SI[, 4] <- if (sum(lambda[4:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:3]), lambda[4] / sum(lambda[4:20])) else 0
n_SI[, 5] <- if (sum(lambda[5:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:4]), lambda[5] / sum(lambda[5:20])) else 0
n_SI[, 6] <- if (sum(lambda[6:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:5]), lambda[6] / sum(lambda[6:20])) else 0
n_SI[, 7] <- if (sum(lambda[7:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:6]), lambda[7] / sum(lambda[7:20])) else 0
n_SI[, 8] <- if (sum(lambda[8:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:7]), lambda[8] / sum(lambda[8:20])) else 0
n_SI[, 9] <- if (sum(lambda[9:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:8]), lambda[9] / sum(lambda[9:20])) else 0
n_SI[, 10] <- if (sum(lambda[10:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:9]), lambda[10] / sum(lambda[10:20])) else 0
n_SI[, 11] <- if (sum(lambda[11:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:10]), lambda[11] / sum(lambda[11:20])) else 0
n_SI[, 12] <- if (sum(lambda[12:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:11]), lambda[12] / sum(lambda[12:20])) else 0
n_SI[, 13] <- if (sum(lambda[13:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:12]), lambda[13] / sum(lambda[13:20])) else 0
n_SI[, 14] <- if (sum(lambda[14:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:13]), lambda[14] / sum(lambda[14:20])) else 0
n_SI[, 15] <- if (sum(lambda[15:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:14]), lambda[15] / sum(lambda[15:20])) else 0
n_SI[, 16] <- if (sum(lambda[16:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:15]), lambda[16] / sum(lambda[16:20])) else 0
n_SI[, 17] <- if (sum(lambda[17:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:16]), lambda[17] / sum(lambda[17:20])) else 0
n_SI[, 18] <- if (sum(lambda[18:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:17]), lambda[18] / sum(lambda[18:20])) else 0
n_SI[, 19] <- if (sum(lambda[19:20]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:18]), lambda[19] / sum(lambda[19:20])) else 0
n_SI[, 20] <- S_out[i] - sum(n_SI[i, 1:19])

## After 1 infection 
### Can be infected with 2,3 & 4
n_RI[, 1, 6] <- if ((sum(lambda[6:20])) > 0) Binomial(R_out[i, 1], lambda[6] / (sum(lambda[6:20]))) else 0
n_RI[, 1, 7] <- if ((sum(lambda[7:20])) > 0) Binomial(R_out[i, 1] - (n_RI[i, 1, 6]), lambda[7] / (sum(lambda[7:20]))) else 0
n_RI[, 1, 8] <- if ((sum(lambda[8:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:7])), lambda[8] / (sum(lambda[8:20]))) else 0
n_RI[, 1, 9] <- if ((sum(lambda[9:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:8])), lambda[9] / (sum(lambda[9:20]))) else 0
n_RI[, 1, 10] <- if ((sum(lambda[10:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:9])), lambda[10] / (sum(lambda[10:20]))) else 0
n_RI[, 1, 11] <- if ((sum(lambda[11:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:10])), lambda[11] / (sum(lambda[11:20]))) else 0
n_RI[, 1, 12] <- if ((sum(lambda[12:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:11])), lambda[12] / (sum(lambda[12:20]))) else 0
n_RI[, 1, 13] <- if ((sum(lambda[13:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:12])), lambda[13] / (sum(lambda[13:20]))) else 0
n_RI[, 1, 14] <- if ((sum(lambda[14:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:13])), lambda[14] / (sum(lambda[14:20]))) else 0
n_RI[, 1, 15] <- if ((sum(lambda[15:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:14])), lambda[15] / (sum(lambda[15:20]))) else 0
n_RI[, 1, 16] <- if ((sum(lambda[16:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:15])), lambda[16] / (sum(lambda[16:20]))) else 0
n_RI[, 1, 17] <- if ((sum(lambda[17:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:16])), lambda[17] / (sum(lambda[17:20]))) else 0
n_RI[, 1, 18] <- if ((sum(lambda[18:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:17])), lambda[18] / (sum(lambda[18:20]))) else 0
n_RI[, 1, 19] <- if ((sum(lambda[19:20])) > 0) Binomial(R_out[i, 1] - (sum(n_RI[i, 1, 6:18])), lambda[19] / (sum(lambda[19:20]))) else 0
n_RI[, 1, 20] <- R_out[i, 1] - (sum(n_RI[i, 1, 6:19]))
n_RI[, 2, 1] <- if ((sum(lambda[1:5]) + sum(lambda[11:20])) > 0) Binomial(R_out[i, 2], lambda[1] / (sum(lambda[1:5]) + sum(lambda[11:20]))) else 0
n_RI[, 2, 2] <- if ((sum(lambda[2:5]) + sum(lambda[11:20])) > 0) Binomial(R_out[i, 2] - (n_RI[i, 2, 1]), lambda[2] / (sum(lambda[2:5]) + sum(lambda[11:20]))) else 0
n_RI[, 2, 3] <- if ((sum(lambda[3:5]) + sum(lambda[11:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:2])), lambda[3] / (sum(lambda[3:5]) + sum(lambda[11:20]))) else 0
n_RI[, 2, 4] <- if ((sum(lambda[4:5]) + sum(lambda[11:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:3])), lambda[4] / (sum(lambda[4:5]) + sum(lambda[11:20]))) else 0
n_RI[, 2, 5] <- if ((lambda[5] + sum(lambda[11:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:4])), lambda[5] / (lambda[5] + sum(lambda[11:20]))) else 0
n_RI[, 2, 11] <- if ((sum(lambda[11:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:5])), lambda[11] / (sum(lambda[11:20]))) else 0
n_RI[, 2, 12] <- if ((sum(lambda[12:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:5]) + n_RI[i, 2, 11]), lambda[12] / (sum(lambda[12:20]))) else 0
n_RI[, 2, 13] <- if ((sum(lambda[13:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:5]) + sum(n_RI[i, 2, 11:12])), lambda[13] / (sum(lambda[13:20]))) else 0
n_RI[, 2, 14] <- if ((sum(lambda[14:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:5]) + sum(n_RI[i, 2, 11:13])), lambda[14] / (sum(lambda[14:20]))) else 0
n_RI[, 2, 15] <- if ((sum(lambda[15:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:5]) + sum(n_RI[i, 2, 11:14])), lambda[15] / (sum(lambda[15:20]))) else 0
n_RI[, 2, 16] <- if ((sum(lambda[16:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:5]) + sum(n_RI[i, 2, 11:15])), lambda[16] / (sum(lambda[16:20]))) else 0
n_RI[, 2, 17] <- if ((sum(lambda[17:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:5]) + sum(n_RI[i, 2, 11:16])), lambda[17] / (sum(lambda[17:20]))) else 0
n_RI[, 2, 18] <- if ((sum(lambda[18:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:5]) + sum(n_RI[i, 2, 11:17])), lambda[18] / (sum(lambda[18:20]))) else 0
n_RI[, 2, 19] <- if ((sum(lambda[19:20])) > 0) Binomial(R_out[i, 2] - (sum(n_RI[i, 2, 1:5]) + sum(n_RI[i, 2, 11:18])), lambda[19] / (sum(lambda[19:20]))) else 0
n_RI[, 2, 20] <- R_out[i, 2] - (sum(n_RI[i, 2, 1:5]) + sum(n_RI[i, 2, 11:19]))
n_RI[, 3, 1] <- if ((sum(lambda[1:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 3], lambda[1] / (sum(lambda[1:10]) + sum(lambda[16:20]))) else 0
n_RI[, 3, 2] <- if ((sum(lambda[2:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 3] - (n_RI[i, 3, 1]), lambda[2] / (sum(lambda[2:10]) + sum(lambda[16:20]))) else 0
n_RI[, 3, 3] <- if ((sum(lambda[3:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:2])), lambda[3] / (sum(lambda[3:10]) + sum(lambda[16:20]))) else 0
n_RI[, 3, 4] <- if ((sum(lambda[4:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:3])), lambda[4] / (sum(lambda[4:10]) + sum(lambda[16:20]))) else 0
n_RI[, 3, 5] <- if ((sum(lambda[5:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:4])), lambda[5] / (sum(lambda[5:10]) + sum(lambda[16:20]))) else 0
n_RI[, 3, 6] <- if ((sum(lambda[6:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:5])), lambda[6] / (sum(lambda[6:10]) + sum(lambda[16:20]))) else 0
n_RI[, 3, 7] <- if ((sum(lambda[7:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:6])), lambda[7] / (sum(lambda[7:10]) + sum(lambda[16:20]))) else 0
n_RI[, 3, 8] <- if ((sum(lambda[8:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:7])), lambda[8] / (sum(lambda[8:10]) + sum(lambda[16:20]))) else 0
n_RI[, 3, 9] <- if ((sum(lambda[9:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:8])), lambda[9] / (sum(lambda[9:10]) + sum(lambda[16:20]))) else 0
n_RI[, 3, 10] <- if ((lambda[10] + sum(lambda[16:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:9])), lambda[10] / (lambda[10] + sum(lambda[16:20]))) else 0
n_RI[, 3, 16] <- if ((sum(lambda[16:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:10])), lambda[16] / (sum(lambda[16:20]))) else 0
n_RI[, 3, 17] <- if ((sum(lambda[17:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:10]) + n_RI[i, 3, 16]), lambda[17] / (sum(lambda[17:20]))) else 0
n_RI[, 3, 18] <- if ((sum(lambda[18:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:10]) + sum(n_RI[i, 3, 16:17])), lambda[18] / (sum(lambda[18:20]))) else 0
n_RI[, 3, 19] <- if ((sum(lambda[19:20])) > 0) Binomial(R_out[i, 3] - (sum(n_RI[i, 3, 1:10]) + sum(n_RI[i, 3, 16:18])), lambda[19] / (sum(lambda[19:20]))) else 0
n_RI[, 3, 20] <- R_out[i, 3] - (sum(n_RI[i, 3, 1:10]) + sum(n_RI[i, 3, 16:19]))
n_RI[, 4, 1] <- if ((sum(lambda[1:15])) > 0) Binomial(R_out[i, 4], lambda[1] / (sum(lambda[1:15]))) else 0
n_RI[, 4, 2] <- if ((sum(lambda[2:15])) > 0) Binomial(R_out[i, 4] - (n_RI[i, 4, 1]), lambda[2] / (sum(lambda[2:15]))) else 0
n_RI[, 4, 3] <- if ((sum(lambda[3:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:2])), lambda[3] / (sum(lambda[3:15]))) else 0
n_RI[, 4, 4] <- if ((sum(lambda[4:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:3])), lambda[4] / (sum(lambda[4:15]))) else 0
n_RI[, 4, 5] <- if ((sum(lambda[5:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:4])), lambda[5] / (sum(lambda[5:15]))) else 0
n_RI[, 4, 6] <- if ((sum(lambda[6:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:5])), lambda[6] / (sum(lambda[6:15]))) else 0
n_RI[, 4, 7] <- if ((sum(lambda[7:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:6])), lambda[7] / (sum(lambda[7:15]))) else 0
n_RI[, 4, 8] <- if ((sum(lambda[8:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:7])), lambda[8] / (sum(lambda[8:15]))) else 0
n_RI[, 4, 9] <- if ((sum(lambda[9:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:8])), lambda[9] / (sum(lambda[9:15]))) else 0
n_RI[, 4, 10] <- if ((sum(lambda[10:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:9])), lambda[10] / (sum(lambda[10:15]))) else 0
n_RI[, 4, 11] <- if ((sum(lambda[11:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:10])), lambda[11] / (sum(lambda[11:15]))) else 0
n_RI[, 4, 12] <- if ((sum(lambda[12:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:11])), lambda[12] / (sum(lambda[12:15]))) else 0
n_RI[, 4, 13] <- if ((sum(lambda[13:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:12])), lambda[13] / (sum(lambda[13:15]))) else 0
n_RI[, 4, 14] <- if ((sum(lambda[14:15])) > 0) Binomial(R_out[i, 4] - (sum(n_RI[i, 4, 1:13])), lambda[14] / (sum(lambda[14:15]))) else 0
n_RI[, 4, 15] <- R_out[i, 4] - (sum(n_RI[i, 4, 1:14]))
n_RI[, 5, 11] <- if ((sum(lambda[11:20])) > 0) Binomial(R_out[i, 5], lambda[11] / (sum(lambda[11:20]))) else 0
n_RI[, 5, 12] <- if ((sum(lambda[12:20])) > 0) Binomial(R_out[i, 5] - (n_RI[i, 5, 11]), lambda[12] / (sum(lambda[12:20]))) else 0
n_RI[, 5, 13] <- if ((sum(lambda[13:20])) > 0) Binomial(R_out[i, 5] - (sum(n_RI[i, 5, 11:12])), lambda[13] / (sum(lambda[13:20]))) else 0
n_RI[, 5, 14] <- if ((sum(lambda[14:20])) > 0) Binomial(R_out[i, 5] - (sum(n_RI[i, 5, 11:13])), lambda[14] / (sum(lambda[14:20]))) else 0
n_RI[, 5, 15] <- if ((sum(lambda[15:20])) > 0) Binomial(R_out[i, 5] - (sum(n_RI[i, 5, 11:14])), lambda[15] / (sum(lambda[15:20]))) else 0
n_RI[, 5, 16] <- if ((sum(lambda[16:20])) > 0) Binomial(R_out[i, 5] - (sum(n_RI[i, 5, 11:15])), lambda[16] / (sum(lambda[16:20]))) else 0
n_RI[, 5, 17] <- if ((sum(lambda[17:20])) > 0) Binomial(R_out[i, 5] - (sum(n_RI[i, 5, 11:16])), lambda[17] / (sum(lambda[17:20]))) else 0
n_RI[, 5, 18] <- if ((sum(lambda[18:20])) > 0) Binomial(R_out[i, 5] - (sum(n_RI[i, 5, 11:17])), lambda[18] / (sum(lambda[18:20]))) else 0
n_RI[, 5, 19] <- if ((sum(lambda[19:20])) > 0) Binomial(R_out[i, 5] - (sum(n_RI[i, 5, 11:18])), lambda[19] / (sum(lambda[19:20]))) else 0
n_RI[, 5, 20] <- R_out[i, 5] - (sum(n_RI[i, 5, 11:19]))
n_RI[, 6, 6] <- if ((sum(lambda[6:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 6], lambda[6] / (sum(lambda[6:10]) + sum(lambda[16:20]))) else 0
n_RI[, 6, 7] <- if ((sum(lambda[7:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 6] - (n_RI[i, 6, 6]), lambda[7] / (sum(lambda[7:10]) + sum(lambda[16:20]))) else 0
n_RI[, 6, 8] <- if ((sum(lambda[8:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 6] - (sum(n_RI[i, 6, 6:7])), lambda[8] / (sum(lambda[8:10]) + sum(lambda[16:20]))) else 0
n_RI[, 6, 9] <- if ((sum(lambda[9:10]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 6] - (sum(n_RI[i, 6, 6:8])), lambda[9] / (sum(lambda[9:10]) + sum(lambda[16:20]))) else 0
n_RI[, 6, 10] <- if ((lambda[10] + sum(lambda[16:20])) > 0) Binomial(R_out[i, 6] - (sum(n_RI[i, 6, 6:9])), lambda[10] / (lambda[10] + sum(lambda[16:20]))) else 0
n_RI[, 6, 16] <- if ((sum(lambda[16:20])) > 0) Binomial(R_out[i, 6] - (sum(n_RI[i, 6, 6:10])), lambda[16] / (sum(lambda[16:20]))) else 0
n_RI[, 6, 17] <- if ((sum(lambda[17:20])) > 0) Binomial(R_out[i, 6] - (sum(n_RI[i, 6, 6:10]) + n_RI[i, 6, 16]), lambda[17] / (sum(lambda[17:20]))) else 0
n_RI[, 6, 18] <- if ((sum(lambda[18:20])) > 0) Binomial(R_out[i, 6] - (sum(n_RI[i, 6, 6:10]) + sum(n_RI[i, 6, 16:17])), lambda[18] / (sum(lambda[18:20]))) else 0
n_RI[, 6, 19] <- if ((sum(lambda[19:20])) > 0) Binomial(R_out[i, 6] - (sum(n_RI[i, 6, 6:10]) + sum(n_RI[i, 6, 16:18])), lambda[19] / (sum(lambda[19:20]))) else 0
n_RI[, 6, 20] <- R_out[i, 6] - (sum(n_RI[i, 6, 6:10]) + sum(n_RI[i, 6, 16:19]))
n_RI[, 7, 6] <- if ((sum(lambda[6:15])) > 0) Binomial(R_out[i, 7], lambda[6] / (sum(lambda[6:15]))) else 0
n_RI[, 7, 7] <- if ((sum(lambda[7:15])) > 0) Binomial(R_out[i, 7] - (n_RI[i, 7, 6]), lambda[7] / (sum(lambda[7:15]))) else 0
n_RI[, 7, 8] <- if ((sum(lambda[8:15])) > 0) Binomial(R_out[i, 7] - (sum(n_RI[i, 7, 6:7])), lambda[8] / (sum(lambda[8:15]))) else 0
n_RI[, 7, 9] <- if ((sum(lambda[9:15])) > 0) Binomial(R_out[i, 7] - (sum(n_RI[i, 7, 6:8])), lambda[9] / (sum(lambda[9:15]))) else 0
n_RI[, 7, 10] <- if ((sum(lambda[10:15])) > 0) Binomial(R_out[i, 7] - (sum(n_RI[i, 7, 6:9])), lambda[10] / (sum(lambda[10:15]))) else 0
n_RI[, 7, 11] <- if ((sum(lambda[11:15])) > 0) Binomial(R_out[i, 7] - (sum(n_RI[i, 7, 6:10])), lambda[11] / (sum(lambda[11:15]))) else 0
n_RI[, 7, 12] <- if ((sum(lambda[12:15])) > 0) Binomial(R_out[i, 7] - (sum(n_RI[i, 7, 6:11])), lambda[12] / (sum(lambda[12:15]))) else 0
n_RI[, 7, 13] <- if ((sum(lambda[13:15])) > 0) Binomial(R_out[i, 7] - (sum(n_RI[i, 7, 6:12])), lambda[13] / (sum(lambda[13:15]))) else 0
n_RI[, 7, 14] <- if ((sum(lambda[14:15])) > 0) Binomial(R_out[i, 7] - (sum(n_RI[i, 7, 6:13])), lambda[14] / (sum(lambda[14:15]))) else 0
n_RI[, 7, 15] <- R_out[i, 7] - (sum(n_RI[i, 7, 6:14]))
n_RI[, 8, 1] <- if ((sum(lambda[1:5]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 8], lambda[1] / (sum(lambda[1:5]) + sum(lambda[16:20]))) else 0
n_RI[, 8, 2] <- if ((sum(lambda[2:5]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 8] - (n_RI[i, 8, 1]), lambda[2] / (sum(lambda[2:5]) + sum(lambda[16:20]))) else 0
n_RI[, 8, 3] <- if ((sum(lambda[3:5]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 8] - (sum(n_RI[i, 8, 1:2])), lambda[3] / (sum(lambda[3:5]) + sum(lambda[16:20]))) else 0
n_RI[, 8, 4] <- if ((sum(lambda[4:5]) + sum(lambda[16:20])) > 0) Binomial(R_out[i, 8] - (sum(n_RI[i, 8, 1:3])), lambda[4] / (sum(lambda[4:5]) + sum(lambda[16:20]))) else 0
n_RI[, 8, 5] <- if ((lambda[5] + sum(lambda[16:20])) > 0) Binomial(R_out[i, 8] - (sum(n_RI[i, 8, 1:4])), lambda[5] / (lambda[5] + sum(lambda[16:20]))) else 0
n_RI[, 8, 16] <- if ((sum(lambda[16:20])) > 0) Binomial(R_out[i, 8] - (sum(n_RI[i, 8, 1:5])), lambda[16] / (sum(lambda[16:20]))) else 0
n_RI[, 8, 17] <- if ((sum(lambda[17:20])) > 0) Binomial(R_out[i, 8] - (sum(n_RI[i, 8, 1:5]) + n_RI[i, 8, 16]), lambda[17] / (sum(lambda[17:20]))) else 0
n_RI[, 8, 18] <- if ((sum(lambda[18:20])) > 0) Binomial(R_out[i, 8] - (sum(n_RI[i, 8, 1:5]) + sum(n_RI[i, 8, 16:17])), lambda[18] / (sum(lambda[18:20]))) else 0
n_RI[, 8, 19] <- if ((sum(lambda[19:20])) > 0) Binomial(R_out[i, 8] - (sum(n_RI[i, 8, 1:5]) + sum(n_RI[i, 8, 16:18])), lambda[19] / (sum(lambda[19:20]))) else 0
n_RI[, 8, 20] <- R_out[i, 8] - (sum(n_RI[i, 8, 1:5]) + sum(n_RI[i, 8, 16:19]))
n_RI[, 9, 1] <- if ((sum(lambda[1:5]) + sum(lambda[11:15])) > 0) Binomial(R_out[i, 9], lambda[1] / (sum(lambda[1:5]) + sum(lambda[11:15]))) else 0
n_RI[, 9, 2] <- if ((sum(lambda[2:5]) + sum(lambda[11:15])) > 0) Binomial(R_out[i, 9] - (n_RI[i, 9, 1]), lambda[2] / (sum(lambda[2:5]) + sum(lambda[11:15]))) else 0
n_RI[, 9, 3] <- if ((sum(lambda[3:5]) + sum(lambda[11:15])) > 0) Binomial(R_out[i, 9] - (sum(n_RI[i, 9, 1:2])), lambda[3] / (sum(lambda[3:5]) + sum(lambda[11:15]))) else 0
n_RI[, 9, 4] <- if ((sum(lambda[4:5]) + sum(lambda[11:15])) > 0) Binomial(R_out[i, 9] - (sum(n_RI[i, 9, 1:3])), lambda[4] / (sum(lambda[4:5]) + sum(lambda[11:15]))) else 0
n_RI[, 9, 5] <- if ((lambda[5] + sum(lambda[11:15])) > 0) Binomial(R_out[i, 9] - (sum(n_RI[i, 9, 1:4])), lambda[5] / (lambda[5] + sum(lambda[11:15]))) else 0
n_RI[, 9, 11] <- if ((sum(lambda[11:15])) > 0) Binomial(R_out[i, 9] - (sum(n_RI[i, 9, 1:5])), lambda[11] / (sum(lambda[11:15]))) else 0
n_RI[, 9, 12] <- if ((sum(lambda[12:15])) > 0) Binomial(R_out[i, 9] - (sum(n_RI[i, 9, 1:5]) + n_RI[i, 9, 11]), lambda[12] / (sum(lambda[12:15]))) else 0
n_RI[, 9, 13] <- if ((sum(lambda[13:15])) > 0) Binomial(R_out[i, 9] - (sum(n_RI[i, 9, 1:5]) + sum(n_RI[i, 9, 11:12])), lambda[13] / (sum(lambda[13:15]))) else 0
n_RI[, 9, 14] <- if ((sum(lambda[14:15])) > 0) Binomial(R_out[i, 9] - (sum(n_RI[i, 9, 1:5]) + sum(n_RI[i, 9, 11:13])), lambda[14] / (sum(lambda[14:15]))) else 0
n_RI[, 9, 15] <- R_out[i, 9] - (sum(n_RI[i, 9, 1:5]) + sum(n_RI[i, 9, 11:14]))
n_RI[, 10, 1] <- if ((sum(lambda[1:10])) > 0) Binomial(R_out[i, 10], lambda[1] / (sum(lambda[1:10]))) else 0
n_RI[, 10, 2] <- if ((sum(lambda[2:10])) > 0) Binomial(R_out[i, 10] - (n_RI[i, 10, 1]), lambda[2] / (sum(lambda[2:10]))) else 0
n_RI[, 10, 3] <- if ((sum(lambda[3:10])) > 0) Binomial(R_out[i, 10] - (sum(n_RI[i, 10, 1:2])), lambda[3] / (sum(lambda[3:10]))) else 0
n_RI[, 10, 4] <- if ((sum(lambda[4:10])) > 0) Binomial(R_out[i, 10] - (sum(n_RI[i, 10, 1:3])), lambda[4] / (sum(lambda[4:10]))) else 0
n_RI[, 10, 5] <- if ((sum(lambda[5:10])) > 0) Binomial(R_out[i, 10] - (sum(n_RI[i, 10, 1:4])), lambda[5] / (sum(lambda[5:10]))) else 0
n_RI[, 10, 6] <- if ((sum(lambda[6:10])) > 0) Binomial(R_out[i, 10] - (sum(n_RI[i, 10, 1:5])), lambda[6] / (sum(lambda[6:10]))) else 0
n_RI[, 10, 7] <- if ((sum(lambda[7:10])) > 0) Binomial(R_out[i, 10] - (sum(n_RI[i, 10, 1:6])), lambda[7] / (sum(lambda[7:10]))) else 0
n_RI[, 10, 8] <- if ((sum(lambda[8:10])) > 0) Binomial(R_out[i, 10] - (sum(n_RI[i, 10, 1:7])), lambda[8] / (sum(lambda[8:10]))) else 0
n_RI[, 10, 9] <- if ((sum(lambda[9:10])) > 0) Binomial(R_out[i, 10] - (sum(n_RI[i, 10, 1:8])), lambda[9] / (sum(lambda[9:10]))) else 0
n_RI[, 10, 10] <- R_out[i, 10] - (sum(n_RI[i, 10, 1:9]))
n_RI[, 11, 16] <- if ((sum(lambda[16:20])) > 0) Binomial(R_out[i, 11], lambda[16] / (sum(lambda[16:20]))) else 0
n_RI[, 11, 17] <- if ((sum(lambda[17:20])) > 0) Binomial(R_out[i, 11] - (n_RI[i, 11, 16]), lambda[17] / (sum(lambda[17:20]))) else 0
n_RI[, 11, 18] <- if ((sum(lambda[18:20])) > 0) Binomial(R_out[i, 11] - (sum(n_RI[i, 11, 16:17])), lambda[18] / (sum(lambda[18:20]))) else 0
n_RI[, 11, 19] <- if ((sum(lambda[19:20])) > 0) Binomial(R_out[i, 11] - (sum(n_RI[i, 11, 16:18])), lambda[19] / (sum(lambda[19:20]))) else 0
n_RI[, 11, 20] <- R_out[i, 11] - (sum(n_RI[i, 11, 16:19]))
n_RI[, 12, 11] <- if ((sum(lambda[11:15])) > 0) Binomial(R_out[i, 12], lambda[11] / (sum(lambda[11:15]))) else 0
n_RI[, 12, 12] <- if ((sum(lambda[12:15])) > 0) Binomial(R_out[i, 12] - (n_RI[i, 12, 11]), lambda[12] / (sum(lambda[12:15]))) else 0
n_RI[, 12, 13] <- if ((sum(lambda[13:15])) > 0) Binomial(R_out[i, 12] - (sum(n_RI[i, 12, 11:12])), lambda[13] / (sum(lambda[13:15]))) else 0
n_RI[, 12, 14] <- if ((sum(lambda[14:15])) > 0) Binomial(R_out[i, 12] - (sum(n_RI[i, 12, 11:13])), lambda[14] / (sum(lambda[14:15]))) else 0
n_RI[, 12, 15] <- R_out[i, 12] - (sum(n_RI[i, 12, 11:14]))
n_RI[, 13, 6] <- if ((sum(lambda[6:10])) > 0) Binomial(R_out[i, 13], lambda[6] / (sum(lambda[6:10]))) else 0
n_RI[, 13, 7] <- if ((sum(lambda[7:10])) > 0) Binomial(R_out[i, 13] - (n_RI[i, 13, 6]), lambda[7] / (sum(lambda[7:10]))) else 0
n_RI[, 13, 8] <- if ((sum(lambda[8:10])) > 0) Binomial(R_out[i, 13] - (sum(n_RI[i, 13, 6:7])), lambda[8] / (sum(lambda[8:10]))) else 0
n_RI[, 13, 9] <- if ((sum(lambda[9:10])) > 0) Binomial(R_out[i, 13] - (sum(n_RI[i, 13, 6:8])), lambda[9] / (sum(lambda[9:10]))) else 0
n_RI[, 13, 10] <- R_out[i, 13] - (sum(n_RI[i, 13, 6:9]))
n_RI[, 14, 1] <- if ((sum(lambda[1:5])) > 0) Binomial(R_out[i, 14], lambda[1] / (sum(lambda[1:5]))) else 0
n_RI[, 14, 2] <- if ((sum(lambda[2:5])) > 0) Binomial(R_out[i, 14] - (n_RI[i, 14, 1]), lambda[2] / (sum(lambda[2:5]))) else 0
n_RI[, 14, 3] <- if ((sum(lambda[3:5])) > 0) Binomial(R_out[i, 14] - (sum(n_RI[i, 14, 1:2])), lambda[3] / (sum(lambda[3:5]))) else 0
n_RI[, 14, 4] <- if ((sum(lambda[4:5])) > 0) Binomial(R_out[i, 14] - (sum(n_RI[i, 14, 1:3])), lambda[4] / (sum(lambda[4:5]))) else 0
n_RI[, 14, 5] <- R_out[i, 14] - (sum(n_RI[i, 14, 1:4]))