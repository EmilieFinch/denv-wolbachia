#### Definition of discrete time stochastic compartmental model
### Structure of script
## 1. Model parameters
## 2. Core equations
## 3. Compute probabilities of transition
## 4. Draw numbers moving between compartments
## 5. Ageing & mortality
## 6. Burden tracking
## 7. Initial states and dimensions

### Compartments S, I and R are stratified by 
## i for the age group (number of age groups specified by n_age, with age group width specified in size_age)
## j for state
## Compartments C and R are also stratified by 
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
n_strains <- 60 #total number of strains modelled
r0 <- parameter() # input strain-specific R0 (e.g. each strain has this R0 when introduced alone)
gamma <- parameter(0.2)
N_init <- parameter() # initial population size by age
demog <- parameter() # array of population size [age_group, sim_year]
ageing_day <- (time %% 365 == 0 && time != 0) # population ages once a year
initial(day_of_year, zero_every = 365) <- 0
update(day_of_year) <- day_of_year + 1

## Wolbachia parameters
wol_on <- parameter()
wol_start <- parameter()
wol_inhib <- parameter() # level of wolbachia inhibition (by strain)

#### Core equations ####
update(S[]) <- if(i == 1 && ageing_day) births[day_of_year + 1] else 
  if(i == 1) births[day_of_year] + S[i] - S_out[i] else
  if(ageing_day && i > 1) S[i-1] - S_out[i-1] else
  S[i] - S_out[i]

update(I[,1,]) <- if(ageing_day && i == 1) 0 else
  if(ageing_day && i > 1)  I[i-1,j,k] - I_out[i-1,j,k] + n_SI[i-1,k] else 
    I[i,j,k] - I_out[i,j,k] + n_SI[i,k]

update(I[,2:n_histories,]) <- if(ageing_day && i == 1) 0 else
  if(ageing_day && i > 1) I[i-1,j,k] - I_out[i-1,j,k] + n_RI[i-1,j-1,k] else
    I[i,j,k] - I_out[i,j,k] + n_RI[i,j-1,k]

update(C[,]) <- if(ageing_day && i == 1) 0 else
 if(ageing_day && i > 1) C[i-1,j] - n_CR[i-1,j] + n_IC[i-1,j] else
    C[i,j] - n_CR[i,j] + n_IC[i,j]

update(R[,]) <- if(ageing_day && i == 1) 0 else
  if(ageing_day && i > 1) R[i-1,j] - R_out[i-1,j] + n_CR[i-1,j] else
    R[i,j] - R_out[i,j] + n_CR[i,j]

update(N[]) <- if(ageing_day) demog[i, sim_year + 1] else N[i] # don't use sum of S, I, C and R to avoid introducing lag in N

#### Calculate individual probabilities of transition ####
beta[] <- r0*gamma
beta_adj[] <- if(wol_on == 1 && time >= wol_start) wol_inhib[i]*beta[i] else beta[i]

lambda[] <- beta_adj[i] * sum(I[,,i])/sum(N_init[]) # strain specific lambda
lambda[] <- lambda[i] 
lambda_1 <- sum(lambda[1:20]) # DENV1 lambda
lambda_2 <- sum(lambda[21:40]) # DENV2 lambda
lambda_3 <- sum(lambda[41:50]) # DENV3 lambda
lambda_4 <- sum(lambda[51:60])

lambda_total <- sum(lambda[])

p_IC <- 1 - exp(-gamma*dt) # I to C
p_CR <- 1 - exp(-nu*dt) # C to R
nu <- 1/365  # duration of cross protection of a year 

print("ageing_day: {ageing_day}, day_of_year: {day_of_year}, births: {births[1]}, sim_year: {sim_year}, N:{N[1]}, S[1]:{S[1]}")

#### Draw number moving between compartments ####
## Note chain binomials for n_SI and n_RI at the end of the script 
S_out[] <- Binomial(S[i], 1 - exp(-lambda_total*dt)) 
I_out[,,] <- Binomial(I[i,j,k], p_IC)

n_IC[, 1] <- sum(I_out[i, 1, 1:20]) # seronaive infected with DENV1
n_IC[, 2] <- sum(I_out[i, 1, 21:40]) # seronaive infected with DENV2
n_IC[, 3] <- sum(I_out[i, 1, 41:50]) # seronaive infected with DENV3
n_IC[, 4] <- sum(I_out[i, 1, 51:60]) # seronaive infected with DENV4
n_IC[, 5] <- sum(I_out[i, 2, 21:40]) + sum(I_out[i, 3, 1:20]) # Prev DENV1 infected with DENV2 and vice versa
n_IC[, 6] <- sum(I_out[i, 2, 41:50]) + sum(I_out[i, 4, 1:20]) # Prev DENV1 infected with DENV3 and vice versa
n_IC[, 7] <- sum(I_out[i, 2, 51:60]) + sum(I_out[i, 5, 1:20]) # Prev DENV1 infected with DENV4 and vice versa
n_IC[, 8] <- sum(I_out[i, 3, 41:50]) + sum(I_out[i, 4, 21:40]) # Prev DENV2 infected with DENV3 and vice versa
n_IC[, 9] <- sum(I_out[i, 3, 51:60]) + sum(I_out[i, 5, 21:40]) # Prev DENV2 infected with DENV4 and vice versa
n_IC[, 10] <- sum(I_out[i, 4, 51:60]) + sum(I_out[i, 5, 41:50]) # Prev DENV3 infected with DENV4 and vice versa
n_IC[, 11] <- sum(I_out[i, 6, 41:50]) + sum(I_out[i, 7, 21:40]) + sum(I_out[i, 9, 1:20]) # Prev 2&3 infected with 1, Prev 1&2 infected with 3 and prev 1&3 infected with 2
n_IC[, 12] <- sum(I_out[i, 6, 51:60]) + sum(I_out[i, 8, 21:40]) + sum(I_out[i, 10, 1:20]) # Prev 2&4 infected with 1, Prev 1&4 infected with 2 and prev 1&2 infected with 4
n_IC[, 13] <- sum(I_out[i, 7, 51:60]) + sum(I_out[i, 8, 41:50]) + sum(I_out[i, 11, 1:20]) # Prev 3&4 infected with 1, Prev 1&4 infected with 3 and prev 1&3 infected with 4
n_IC[, 14] <- sum(I_out[i, 9, 51:60]) + sum(I_out[i, 10, 41:50]) + sum(I_out[i, 11, 21:40]) # Prev 3&4 infected with 2, Prev 2&4 infected with 3 and prev 2&3 infected with 4
n_IC[, 15] <- sum(I_out[i, 12, 51:60]) + sum(I_out[i, 13, 41:50]) + sum(I_out[i, 14, 21:40]) + sum(I_out[i, 15, 1:20]) # Prev 1&2&3 infected with 4, prev 1&3&4 infected with 2, prev 2&3&4 infected with 1
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
#current_N[] <- S[i] + sum(I[i,,]) + sum(C[i,]) + sum(R[i,])
#shift[2:n_age] <- if(ageing_day && current_N[i-1] > 0) demog[i, sim_year+1]/current_N[i-1] else 1.0 # adjust population sizes so they match demographic data from upcoming year
##shift[1] <- 1 # for youngest age group all enter as susceptible

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

#### Initial states & dimensions ####
initial(N[]) <- N_init[i]
initial(S[]) <- if(i == 11) N_init[i] - 60 else N_init[i]
initial(I[11,1,1:20]) <- 1
initial(I[11,2,21:40]) <- 1
initial(I[11,3,41:50]) <- 1
initial(I[11,4,51:60]) <- 1

initial(C[,]) <- 0
initial(R[,]) <- 0
initial(inf[]) <- 0
initial(inf_age[]) <- 0
initial(prior_infection[]) <- 0

dim(S) <- n_age
dim(I) <- c(n_age, n_histories, n_strains)
dim(C) <- c(n_age, n_histories)
dim(R) <- c(n_age, n_histories)
dim(S_out) <- n_age
dim(N_init) <- n_age
dim(beta) <- n_strains
dim(beta_adj) <- n_strains
dim(wol_inhib) <- n_strains
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
#dim(shift) <- c(n_age)
#dim(current_N) <- c(n_age)

## Chain binomial for strain specific infections ## 
## From susceptible
n_SI[, 1] <- if (sum(lambda[1:60]) > 0) Binomial(S_out[i], lambda[1] / sum(lambda[1:60])) else 0
n_SI[, 2] <- if (sum(lambda[2:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:1]), lambda[2] / sum(lambda[2:60])) else 0
n_SI[, 3] <- if (sum(lambda[3:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:2]), lambda[3] / sum(lambda[3:60])) else 0
n_SI[, 4] <- if (sum(lambda[4:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:3]), lambda[4] / sum(lambda[4:60])) else 0
n_SI[, 5] <- if (sum(lambda[5:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:4]), lambda[5] / sum(lambda[5:60])) else 0
n_SI[, 6] <- if (sum(lambda[6:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:5]), lambda[6] / sum(lambda[6:60])) else 0
n_SI[, 7] <- if (sum(lambda[7:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:6]), lambda[7] / sum(lambda[7:60])) else 0
n_SI[, 8] <- if (sum(lambda[8:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:7]), lambda[8] / sum(lambda[8:60])) else 0
n_SI[, 9] <- if (sum(lambda[9:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:8]), lambda[9] / sum(lambda[9:60])) else 0
n_SI[, 10] <- if (sum(lambda[10:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:9]), lambda[10] / sum(lambda[10:60])) else 0
n_SI[, 11] <- if (sum(lambda[11:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:10]), lambda[11] / sum(lambda[11:60])) else 0
n_SI[, 12] <- if (sum(lambda[12:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:11]), lambda[12] / sum(lambda[12:60])) else 0
n_SI[, 13] <- if (sum(lambda[13:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:12]), lambda[13] / sum(lambda[13:60])) else 0
n_SI[, 14] <- if (sum(lambda[14:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:13]), lambda[14] / sum(lambda[14:60])) else 0
n_SI[, 15] <- if (sum(lambda[15:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:14]), lambda[15] / sum(lambda[15:60])) else 0
n_SI[, 16] <- if (sum(lambda[16:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:15]), lambda[16] / sum(lambda[16:60])) else 0
n_SI[, 17] <- if (sum(lambda[17:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:16]), lambda[17] / sum(lambda[17:60])) else 0
n_SI[, 18] <- if (sum(lambda[18:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:17]), lambda[18] / sum(lambda[18:60])) else 0
n_SI[, 19] <- if (sum(lambda[19:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:18]), lambda[19] / sum(lambda[19:60])) else 0
n_SI[, 20] <- if (sum(lambda[20:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:19]), lambda[20] / sum(lambda[20:60])) else 0
n_SI[, 21] <- if (sum(lambda[21:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:20]), lambda[21] / sum(lambda[21:60])) else 0
n_SI[, 22] <- if (sum(lambda[22:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:21]), lambda[22] / sum(lambda[22:60])) else 0
n_SI[, 23] <- if (sum(lambda[23:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:22]), lambda[23] / sum(lambda[23:60])) else 0
n_SI[, 24] <- if (sum(lambda[24:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:23]), lambda[24] / sum(lambda[24:60])) else 0
n_SI[, 25] <- if (sum(lambda[25:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:24]), lambda[25] / sum(lambda[25:60])) else 0
n_SI[, 26] <- if (sum(lambda[26:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:25]), lambda[26] / sum(lambda[26:60])) else 0
n_SI[, 27] <- if (sum(lambda[27:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:26]), lambda[27] / sum(lambda[27:60])) else 0
n_SI[, 28] <- if (sum(lambda[28:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:27]), lambda[28] / sum(lambda[28:60])) else 0
n_SI[, 29] <- if (sum(lambda[29:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:28]), lambda[29] / sum(lambda[29:60])) else 0
n_SI[, 30] <- if (sum(lambda[30:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:29]), lambda[30] / sum(lambda[30:60])) else 0
n_SI[, 31] <- if (sum(lambda[31:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:30]), lambda[31] / sum(lambda[31:60])) else 0
n_SI[, 32] <- if (sum(lambda[32:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:31]), lambda[32] / sum(lambda[32:60])) else 0
n_SI[, 33] <- if (sum(lambda[33:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:32]), lambda[33] / sum(lambda[33:60])) else 0
n_SI[, 34] <- if (sum(lambda[34:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:33]), lambda[34] / sum(lambda[34:60])) else 0
n_SI[, 35] <- if (sum(lambda[35:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:34]), lambda[35] / sum(lambda[35:60])) else 0
n_SI[, 36] <- if (sum(lambda[36:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:35]), lambda[36] / sum(lambda[36:60])) else 0
n_SI[, 37] <- if (sum(lambda[37:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:36]), lambda[37] / sum(lambda[37:60])) else 0
n_SI[, 38] <- if (sum(lambda[38:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:37]), lambda[38] / sum(lambda[38:60])) else 0
n_SI[, 39] <- if (sum(lambda[39:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:38]), lambda[39] / sum(lambda[39:60])) else 0
n_SI[, 40] <- if (sum(lambda[40:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:39]), lambda[40] / sum(lambda[40:60])) else 0
n_SI[, 41] <- if (sum(lambda[41:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:40]), lambda[41] / sum(lambda[41:60])) else 0
n_SI[, 42] <- if (sum(lambda[42:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:41]), lambda[42] / sum(lambda[42:60])) else 0
n_SI[, 43] <- if (sum(lambda[43:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:42]), lambda[43] / sum(lambda[43:60])) else 0
n_SI[, 44] <- if (sum(lambda[44:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:43]), lambda[44] / sum(lambda[44:60])) else 0
n_SI[, 45] <- if (sum(lambda[45:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:44]), lambda[45] / sum(lambda[45:60])) else 0
n_SI[, 46] <- if (sum(lambda[46:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:45]), lambda[46] / sum(lambda[46:60])) else 0
n_SI[, 47] <- if (sum(lambda[47:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:46]), lambda[47] / sum(lambda[47:60])) else 0
n_SI[, 48] <- if (sum(lambda[48:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:47]), lambda[48] / sum(lambda[48:60])) else 0
n_SI[, 49] <- if (sum(lambda[49:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:48]), lambda[49] / sum(lambda[49:60])) else 0
n_SI[, 50] <- if (sum(lambda[50:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:49]), lambda[50] / sum(lambda[50:60])) else 0
n_SI[, 51] <- if (sum(lambda[51:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:50]), lambda[51] / sum(lambda[51:60])) else 0
n_SI[, 52] <- if (sum(lambda[52:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:51]), lambda[52] / sum(lambda[52:60])) else 0
n_SI[, 53] <- if (sum(lambda[53:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:52]), lambda[53] / sum(lambda[53:60])) else 0
n_SI[, 54] <- if (sum(lambda[54:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:53]), lambda[54] / sum(lambda[54:60])) else 0
n_SI[, 55] <- if (sum(lambda[55:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:54]), lambda[55] / sum(lambda[55:60])) else 0
n_SI[, 56] <- if (sum(lambda[56:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:55]), lambda[56] / sum(lambda[56:60])) else 0
n_SI[, 57] <- if (sum(lambda[57:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:56]), lambda[57] / sum(lambda[57:60])) else 0
n_SI[, 58] <- if (sum(lambda[58:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:57]), lambda[58] / sum(lambda[58:60])) else 0
n_SI[, 59] <- if (sum(lambda[59:60]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:58]), lambda[59] / sum(lambda[59:60])) else 0
n_SI[, 60] <- S_out[i] - sum(n_SI[i, 1:59])

## After 1 infection

n_RI[, 1, 21] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1], lambda[21] / sum(lambda[21:60])) else 0
n_RI[, 1, 22] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:21]), lambda[22] / sum(lambda[21:60])) else 0
n_RI[, 1, 23] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:22]), lambda[23] / sum(lambda[21:60])) else 0
n_RI[, 1, 24] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:23]), lambda[24] / sum(lambda[21:60])) else 0
n_RI[, 1, 25] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:24]), lambda[25] / sum(lambda[21:60])) else 0
n_RI[, 1, 26] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:25]), lambda[26] / sum(lambda[21:60])) else 0
n_RI[, 1, 27] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:26]), lambda[27] / sum(lambda[21:60])) else 0
n_RI[, 1, 28] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:27]), lambda[28] / sum(lambda[21:60])) else 0
n_RI[, 1, 29] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:28]), lambda[29] / sum(lambda[21:60])) else 0
n_RI[, 1, 30] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:29]), lambda[30] / sum(lambda[21:60])) else 0
n_RI[, 1, 31] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:30]), lambda[31] / sum(lambda[21:60])) else 0
n_RI[, 1, 32] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:31]), lambda[32] / sum(lambda[21:60])) else 0
n_RI[, 1, 33] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:32]), lambda[33] / sum(lambda[21:60])) else 0
n_RI[, 1, 34] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:33]), lambda[34] / sum(lambda[21:60])) else 0
n_RI[, 1, 35] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:34]), lambda[35] / sum(lambda[21:60])) else 0
n_RI[, 1, 36] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:35]), lambda[36] / sum(lambda[21:60])) else 0
n_RI[, 1, 37] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:36]), lambda[37] / sum(lambda[21:60])) else 0
n_RI[, 1, 38] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:37]), lambda[38] / sum(lambda[21:60])) else 0
n_RI[, 1, 39] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:38]), lambda[39] / sum(lambda[21:60])) else 0
n_RI[, 1, 40] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:39]), lambda[40] / sum(lambda[21:60])) else 0
n_RI[, 1, 41] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:40]), lambda[41] / sum(lambda[21:60])) else 0
n_RI[, 1, 42] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:41]), lambda[42] / sum(lambda[21:60])) else 0
n_RI[, 1, 43] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:42]), lambda[43] / sum(lambda[21:60])) else 0
n_RI[, 1, 44] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:43]), lambda[44] / sum(lambda[21:60])) else 0
n_RI[, 1, 45] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:44]), lambda[45] / sum(lambda[21:60])) else 0
n_RI[, 1, 46] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:45]), lambda[46] / sum(lambda[21:60])) else 0
n_RI[, 1, 47] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:46]), lambda[47] / sum(lambda[21:60])) else 0
n_RI[, 1, 48] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:47]), lambda[48] / sum(lambda[21:60])) else 0
n_RI[, 1, 49] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:48]), lambda[49] / sum(lambda[21:60])) else 0
n_RI[, 1, 50] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:49]), lambda[50] / sum(lambda[21:60])) else 0
n_RI[, 1, 51] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:50]), lambda[51] / sum(lambda[21:60])) else 0
n_RI[, 1, 52] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:51]), lambda[52] / sum(lambda[21:60])) else 0
n_RI[, 1, 53] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:52]), lambda[53] / sum(lambda[21:60])) else 0
n_RI[, 1, 54] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:53]), lambda[54] / sum(lambda[21:60])) else 0
n_RI[, 1, 55] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:54]), lambda[55] / sum(lambda[21:60])) else 0
n_RI[, 1, 56] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:55]), lambda[56] / sum(lambda[21:60])) else 0
n_RI[, 1, 57] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:56]), lambda[57] / sum(lambda[21:60])) else 0
n_RI[, 1, 58] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:57]), lambda[58] / sum(lambda[21:60])) else 0
n_RI[, 1, 59] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 1] - sum(n_RI[i, 1, 21:58]), lambda[59] / sum(lambda[21:60])) else 0
n_RI[, 1, 60] <- R_out[i, 1] - sum(n_RI[i, 1, 21:59])
n_RI[, 2, 1] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2], lambda[1] / sum(lambda[1:60])) else 0
n_RI[, 2, 2] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:1]), lambda[2] / sum(lambda[1:60])) else 0
n_RI[, 2, 3] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:2]), lambda[3] / sum(lambda[1:60])) else 0
n_RI[, 2, 4] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:3]), lambda[4] / sum(lambda[1:60])) else 0
n_RI[, 2, 5] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:4]), lambda[5] / sum(lambda[1:60])) else 0
n_RI[, 2, 6] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:5]), lambda[6] / sum(lambda[1:60])) else 0
n_RI[, 2, 7] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:6]), lambda[7] / sum(lambda[1:60])) else 0
n_RI[, 2, 8] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:7]), lambda[8] / sum(lambda[1:60])) else 0
n_RI[, 2, 9] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:8]), lambda[9] / sum(lambda[1:60])) else 0
n_RI[, 2, 10] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:9]), lambda[10] / sum(lambda[1:60])) else 0
n_RI[, 2, 11] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:10]), lambda[11] / sum(lambda[1:60])) else 0
n_RI[, 2, 12] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:11]), lambda[12] / sum(lambda[1:60])) else 0
n_RI[, 2, 13] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:12]), lambda[13] / sum(lambda[1:60])) else 0
n_RI[, 2, 14] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:13]), lambda[14] / sum(lambda[1:60])) else 0
n_RI[, 2, 15] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:14]), lambda[15] / sum(lambda[1:60])) else 0
n_RI[, 2, 16] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:15]), lambda[16] / sum(lambda[1:60])) else 0
n_RI[, 2, 17] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:16]), lambda[17] / sum(lambda[1:60])) else 0
n_RI[, 2, 18] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:17]), lambda[18] / sum(lambda[1:60])) else 0
n_RI[, 2, 19] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:18]), lambda[19] / sum(lambda[1:60])) else 0
n_RI[, 2, 20] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:19]), lambda[20] / sum(lambda[1:60])) else 0
n_RI[, 2, 41] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:20]), lambda[41] / sum(lambda[1:60])) else 0
n_RI[, 2, 42] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:41]), lambda[42] / sum(lambda[1:60])) else 0
n_RI[, 2, 43] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:42]), lambda[43] / sum(lambda[1:60])) else 0
n_RI[, 2, 44] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:43]), lambda[44] / sum(lambda[1:60])) else 0
n_RI[, 2, 45] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:44]), lambda[45] / sum(lambda[1:60])) else 0
n_RI[, 2, 46] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:45]), lambda[46] / sum(lambda[1:60])) else 0
n_RI[, 2, 47] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:46]), lambda[47] / sum(lambda[1:60])) else 0
n_RI[, 2, 48] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:47]), lambda[48] / sum(lambda[1:60])) else 0
n_RI[, 2, 49] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:48]), lambda[49] / sum(lambda[1:60])) else 0
n_RI[, 2, 50] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:49]), lambda[50] / sum(lambda[1:60])) else 0
n_RI[, 2, 51] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:50]), lambda[51] / sum(lambda[1:60])) else 0
n_RI[, 2, 52] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:51]), lambda[52] / sum(lambda[1:60])) else 0
n_RI[, 2, 53] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:52]), lambda[53] / sum(lambda[1:60])) else 0
n_RI[, 2, 54] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:53]), lambda[54] / sum(lambda[1:60])) else 0
n_RI[, 2, 55] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:54]), lambda[55] / sum(lambda[1:60])) else 0
n_RI[, 2, 56] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:55]), lambda[56] / sum(lambda[1:60])) else 0
n_RI[, 2, 57] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:56]), lambda[57] / sum(lambda[1:60])) else 0
n_RI[, 2, 58] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:57]), lambda[58] / sum(lambda[1:60])) else 0
n_RI[, 2, 59] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 2] - sum(n_RI[i, 2, 1:58]), lambda[59] / sum(lambda[1:60])) else 0
n_RI[, 2, 60] <- R_out[i, 2] - sum(n_RI[i, 2, 1:59])
n_RI[, 3, 1] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3], lambda[1] / sum(lambda[1:60])) else 0
n_RI[, 3, 2] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:1]), lambda[2] / sum(lambda[1:60])) else 0
n_RI[, 3, 3] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:2]), lambda[3] / sum(lambda[1:60])) else 0
n_RI[, 3, 4] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:3]), lambda[4] / sum(lambda[1:60])) else 0
n_RI[, 3, 5] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:4]), lambda[5] / sum(lambda[1:60])) else 0
n_RI[, 3, 6] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:5]), lambda[6] / sum(lambda[1:60])) else 0
n_RI[, 3, 7] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:6]), lambda[7] / sum(lambda[1:60])) else 0
n_RI[, 3, 8] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:7]), lambda[8] / sum(lambda[1:60])) else 0
n_RI[, 3, 9] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:8]), lambda[9] / sum(lambda[1:60])) else 0
n_RI[, 3, 10] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:9]), lambda[10] / sum(lambda[1:60])) else 0
n_RI[, 3, 11] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:10]), lambda[11] / sum(lambda[1:60])) else 0
n_RI[, 3, 12] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:11]), lambda[12] / sum(lambda[1:60])) else 0
n_RI[, 3, 13] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:12]), lambda[13] / sum(lambda[1:60])) else 0
n_RI[, 3, 14] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:13]), lambda[14] / sum(lambda[1:60])) else 0
n_RI[, 3, 15] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:14]), lambda[15] / sum(lambda[1:60])) else 0
n_RI[, 3, 16] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:15]), lambda[16] / sum(lambda[1:60])) else 0
n_RI[, 3, 17] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:16]), lambda[17] / sum(lambda[1:60])) else 0
n_RI[, 3, 18] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:17]), lambda[18] / sum(lambda[1:60])) else 0
n_RI[, 3, 19] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:18]), lambda[19] / sum(lambda[1:60])) else 0
n_RI[, 3, 20] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:19]), lambda[20] / sum(lambda[1:60])) else 0
n_RI[, 3, 21] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:20]), lambda[21] / sum(lambda[1:60])) else 0
n_RI[, 3, 22] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:21]), lambda[22] / sum(lambda[1:60])) else 0
n_RI[, 3, 23] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:22]), lambda[23] / sum(lambda[1:60])) else 0
n_RI[, 3, 24] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:23]), lambda[24] / sum(lambda[1:60])) else 0
n_RI[, 3, 25] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:24]), lambda[25] / sum(lambda[1:60])) else 0
n_RI[, 3, 26] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:25]), lambda[26] / sum(lambda[1:60])) else 0
n_RI[, 3, 27] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:26]), lambda[27] / sum(lambda[1:60])) else 0
n_RI[, 3, 28] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:27]), lambda[28] / sum(lambda[1:60])) else 0
n_RI[, 3, 29] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:28]), lambda[29] / sum(lambda[1:60])) else 0
n_RI[, 3, 30] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:29]), lambda[30] / sum(lambda[1:60])) else 0
n_RI[, 3, 31] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:30]), lambda[31] / sum(lambda[1:60])) else 0
n_RI[, 3, 32] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:31]), lambda[32] / sum(lambda[1:60])) else 0
n_RI[, 3, 33] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:32]), lambda[33] / sum(lambda[1:60])) else 0
n_RI[, 3, 34] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:33]), lambda[34] / sum(lambda[1:60])) else 0
n_RI[, 3, 35] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:34]), lambda[35] / sum(lambda[1:60])) else 0
n_RI[, 3, 36] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:35]), lambda[36] / sum(lambda[1:60])) else 0
n_RI[, 3, 37] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:36]), lambda[37] / sum(lambda[1:60])) else 0
n_RI[, 3, 38] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:37]), lambda[38] / sum(lambda[1:60])) else 0
n_RI[, 3, 39] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:38]), lambda[39] / sum(lambda[1:60])) else 0
n_RI[, 3, 40] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:39]), lambda[40] / sum(lambda[1:60])) else 0
n_RI[, 3, 51] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:40]), lambda[51] / sum(lambda[1:60])) else 0
n_RI[, 3, 52] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:51]), lambda[52] / sum(lambda[1:60])) else 0
n_RI[, 3, 53] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:52]), lambda[53] / sum(lambda[1:60])) else 0
n_RI[, 3, 54] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:53]), lambda[54] / sum(lambda[1:60])) else 0
n_RI[, 3, 55] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:54]), lambda[55] / sum(lambda[1:60])) else 0
n_RI[, 3, 56] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:55]), lambda[56] / sum(lambda[1:60])) else 0
n_RI[, 3, 57] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:56]), lambda[57] / sum(lambda[1:60])) else 0
n_RI[, 3, 58] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:57]), lambda[58] / sum(lambda[1:60])) else 0
n_RI[, 3, 59] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 3] - sum(n_RI[i, 3, 1:58]), lambda[59] / sum(lambda[1:60])) else 0
n_RI[, 3, 60] <- R_out[i, 3] - sum(n_RI[i, 3, 1:59])
n_RI[, 4, 1] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4], lambda[1] / sum(lambda[1:50])) else 0
n_RI[, 4, 2] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:1]), lambda[2] / sum(lambda[1:50])) else 0
n_RI[, 4, 3] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:2]), lambda[3] / sum(lambda[1:50])) else 0
n_RI[, 4, 4] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:3]), lambda[4] / sum(lambda[1:50])) else 0
n_RI[, 4, 5] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:4]), lambda[5] / sum(lambda[1:50])) else 0
n_RI[, 4, 6] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:5]), lambda[6] / sum(lambda[1:50])) else 0
n_RI[, 4, 7] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:6]), lambda[7] / sum(lambda[1:50])) else 0
n_RI[, 4, 8] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:7]), lambda[8] / sum(lambda[1:50])) else 0
n_RI[, 4, 9] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:8]), lambda[9] / sum(lambda[1:50])) else 0
n_RI[, 4, 10] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:9]), lambda[10] / sum(lambda[1:50])) else 0
n_RI[, 4, 11] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:10]), lambda[11] / sum(lambda[1:50])) else 0
n_RI[, 4, 12] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:11]), lambda[12] / sum(lambda[1:50])) else 0
n_RI[, 4, 13] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:12]), lambda[13] / sum(lambda[1:50])) else 0
n_RI[, 4, 14] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:13]), lambda[14] / sum(lambda[1:50])) else 0
n_RI[, 4, 15] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:14]), lambda[15] / sum(lambda[1:50])) else 0
n_RI[, 4, 16] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:15]), lambda[16] / sum(lambda[1:50])) else 0
n_RI[, 4, 17] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:16]), lambda[17] / sum(lambda[1:50])) else 0
n_RI[, 4, 18] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:17]), lambda[18] / sum(lambda[1:50])) else 0
n_RI[, 4, 19] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:18]), lambda[19] / sum(lambda[1:50])) else 0
n_RI[, 4, 20] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:19]), lambda[20] / sum(lambda[1:50])) else 0
n_RI[, 4, 21] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:20]), lambda[21] / sum(lambda[1:50])) else 0
n_RI[, 4, 22] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:21]), lambda[22] / sum(lambda[1:50])) else 0
n_RI[, 4, 23] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:22]), lambda[23] / sum(lambda[1:50])) else 0
n_RI[, 4, 24] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:23]), lambda[24] / sum(lambda[1:50])) else 0
n_RI[, 4, 25] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:24]), lambda[25] / sum(lambda[1:50])) else 0
n_RI[, 4, 26] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:25]), lambda[26] / sum(lambda[1:50])) else 0
n_RI[, 4, 27] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:26]), lambda[27] / sum(lambda[1:50])) else 0
n_RI[, 4, 28] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:27]), lambda[28] / sum(lambda[1:50])) else 0
n_RI[, 4, 29] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:28]), lambda[29] / sum(lambda[1:50])) else 0
n_RI[, 4, 30] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:29]), lambda[30] / sum(lambda[1:50])) else 0
n_RI[, 4, 31] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:30]), lambda[31] / sum(lambda[1:50])) else 0
n_RI[, 4, 32] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:31]), lambda[32] / sum(lambda[1:50])) else 0
n_RI[, 4, 33] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:32]), lambda[33] / sum(lambda[1:50])) else 0
n_RI[, 4, 34] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:33]), lambda[34] / sum(lambda[1:50])) else 0
n_RI[, 4, 35] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:34]), lambda[35] / sum(lambda[1:50])) else 0
n_RI[, 4, 36] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:35]), lambda[36] / sum(lambda[1:50])) else 0
n_RI[, 4, 37] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:36]), lambda[37] / sum(lambda[1:50])) else 0
n_RI[, 4, 38] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:37]), lambda[38] / sum(lambda[1:50])) else 0
n_RI[, 4, 39] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:38]), lambda[39] / sum(lambda[1:50])) else 0
n_RI[, 4, 40] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:39]), lambda[40] / sum(lambda[1:50])) else 0
n_RI[, 4, 41] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:40]), lambda[41] / sum(lambda[1:50])) else 0
n_RI[, 4, 42] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:41]), lambda[42] / sum(lambda[1:50])) else 0
n_RI[, 4, 43] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:42]), lambda[43] / sum(lambda[1:50])) else 0
n_RI[, 4, 44] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:43]), lambda[44] / sum(lambda[1:50])) else 0
n_RI[, 4, 45] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:44]), lambda[45] / sum(lambda[1:50])) else 0
n_RI[, 4, 46] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:45]), lambda[46] / sum(lambda[1:50])) else 0
n_RI[, 4, 47] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:46]), lambda[47] / sum(lambda[1:50])) else 0
n_RI[, 4, 48] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:47]), lambda[48] / sum(lambda[1:50])) else 0
n_RI[, 4, 49] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 4] - sum(n_RI[i, 4, 1:48]), lambda[49] / sum(lambda[1:50])) else 0
n_RI[, 4, 50] <- R_out[i, 4] - sum(n_RI[i, 4, 1:49])
n_RI[, 5, 41] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5], lambda[41] / sum(lambda[41:60])) else 0
n_RI[, 5, 42] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:41]), lambda[42] / sum(lambda[41:60])) else 0
n_RI[, 5, 43] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:42]), lambda[43] / sum(lambda[41:60])) else 0
n_RI[, 5, 44] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:43]), lambda[44] / sum(lambda[41:60])) else 0
n_RI[, 5, 45] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:44]), lambda[45] / sum(lambda[41:60])) else 0
n_RI[, 5, 46] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:45]), lambda[46] / sum(lambda[41:60])) else 0
n_RI[, 5, 47] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:46]), lambda[47] / sum(lambda[41:60])) else 0
n_RI[, 5, 48] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:47]), lambda[48] / sum(lambda[41:60])) else 0
n_RI[, 5, 49] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:48]), lambda[49] / sum(lambda[41:60])) else 0
n_RI[, 5, 50] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:49]), lambda[50] / sum(lambda[41:60])) else 0
n_RI[, 5, 51] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:50]), lambda[51] / sum(lambda[41:60])) else 0
n_RI[, 5, 52] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:51]), lambda[52] / sum(lambda[41:60])) else 0
n_RI[, 5, 53] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:52]), lambda[53] / sum(lambda[41:60])) else 0
n_RI[, 5, 54] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:53]), lambda[54] / sum(lambda[41:60])) else 0
n_RI[, 5, 55] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:54]), lambda[55] / sum(lambda[41:60])) else 0
n_RI[, 5, 56] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:55]), lambda[56] / sum(lambda[41:60])) else 0
n_RI[, 5, 57] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:56]), lambda[57] / sum(lambda[41:60])) else 0
n_RI[, 5, 58] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:57]), lambda[58] / sum(lambda[41:60])) else 0
n_RI[, 5, 59] <- if (sum(lambda[41:60]) > 0) Binomial(R_out[i, 5] - sum(n_RI[i, 5, 41:58]), lambda[59] / sum(lambda[41:60])) else 0
n_RI[, 5, 60] <- R_out[i, 5] - sum(n_RI[i, 5, 41:59])

# After 3 infections
n_RI[, 6, 21] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6], lambda[21] / sum(lambda[21:60])) else 0
n_RI[, 6, 22] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:21]), lambda[22] / sum(lambda[21:60])) else 0
n_RI[, 6, 23] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:22]), lambda[23] / sum(lambda[21:60])) else 0
n_RI[, 6, 24] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:23]), lambda[24] / sum(lambda[21:60])) else 0
n_RI[, 6, 25] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:24]), lambda[25] / sum(lambda[21:60])) else 0
n_RI[, 6, 26] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:25]), lambda[26] / sum(lambda[21:60])) else 0
n_RI[, 6, 27] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:26]), lambda[27] / sum(lambda[21:60])) else 0
n_RI[, 6, 28] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:27]), lambda[28] / sum(lambda[21:60])) else 0
n_RI[, 6, 29] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:28]), lambda[29] / sum(lambda[21:60])) else 0
n_RI[, 6, 30] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:29]), lambda[30] / sum(lambda[21:60])) else 0
n_RI[, 6, 31] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:30]), lambda[31] / sum(lambda[21:60])) else 0
n_RI[, 6, 32] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:31]), lambda[32] / sum(lambda[21:60])) else 0
n_RI[, 6, 33] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:32]), lambda[33] / sum(lambda[21:60])) else 0
n_RI[, 6, 34] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:33]), lambda[34] / sum(lambda[21:60])) else 0
n_RI[, 6, 35] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:34]), lambda[35] / sum(lambda[21:60])) else 0
n_RI[, 6, 36] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:35]), lambda[36] / sum(lambda[21:60])) else 0
n_RI[, 6, 37] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:36]), lambda[37] / sum(lambda[21:60])) else 0
n_RI[, 6, 38] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:37]), lambda[38] / sum(lambda[21:60])) else 0
n_RI[, 6, 39] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:38]), lambda[39] / sum(lambda[21:60])) else 0
n_RI[, 6, 40] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:39]), lambda[40] / sum(lambda[21:60])) else 0
n_RI[, 6, 51] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:40]), lambda[51] / sum(lambda[21:60])) else 0
n_RI[, 6, 52] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:51]), lambda[52] / sum(lambda[21:60])) else 0
n_RI[, 6, 53] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:52]), lambda[53] / sum(lambda[21:60])) else 0
n_RI[, 6, 54] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:53]), lambda[54] / sum(lambda[21:60])) else 0
n_RI[, 6, 55] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:54]), lambda[55] / sum(lambda[21:60])) else 0
n_RI[, 6, 56] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:55]), lambda[56] / sum(lambda[21:60])) else 0
n_RI[, 6, 57] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:56]), lambda[57] / sum(lambda[21:60])) else 0
n_RI[, 6, 58] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:57]), lambda[58] / sum(lambda[21:60])) else 0
n_RI[, 6, 59] <- if (sum(lambda[21:60]) > 0) Binomial(R_out[i, 6] - sum(n_RI[i, 6, 21:58]), lambda[59] / sum(lambda[21:60])) else 0
n_RI[, 6, 60] <- R_out[i, 6] - sum(n_RI[i, 6, 21:59])
# After 3 infections
n_RI[, 7, 21] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7], lambda[21] / sum(lambda[21:50])) else 0
n_RI[, 7, 22] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:21]), lambda[22] / sum(lambda[21:50])) else 0
n_RI[, 7, 23] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:22]), lambda[23] / sum(lambda[21:50])) else 0
n_RI[, 7, 24] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:23]), lambda[24] / sum(lambda[21:50])) else 0
n_RI[, 7, 25] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:24]), lambda[25] / sum(lambda[21:50])) else 0
n_RI[, 7, 26] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:25]), lambda[26] / sum(lambda[21:50])) else 0
n_RI[, 7, 27] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:26]), lambda[27] / sum(lambda[21:50])) else 0
n_RI[, 7, 28] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:27]), lambda[28] / sum(lambda[21:50])) else 0
n_RI[, 7, 29] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:28]), lambda[29] / sum(lambda[21:50])) else 0
n_RI[, 7, 30] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:29]), lambda[30] / sum(lambda[21:50])) else 0
n_RI[, 7, 31] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:30]), lambda[31] / sum(lambda[21:50])) else 0
n_RI[, 7, 32] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:31]), lambda[32] / sum(lambda[21:50])) else 0
n_RI[, 7, 33] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:32]), lambda[33] / sum(lambda[21:50])) else 0
n_RI[, 7, 34] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:33]), lambda[34] / sum(lambda[21:50])) else 0
n_RI[, 7, 35] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:34]), lambda[35] / sum(lambda[21:50])) else 0
n_RI[, 7, 36] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:35]), lambda[36] / sum(lambda[21:50])) else 0
n_RI[, 7, 37] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:36]), lambda[37] / sum(lambda[21:50])) else 0
n_RI[, 7, 38] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:37]), lambda[38] / sum(lambda[21:50])) else 0
n_RI[, 7, 39] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:38]), lambda[39] / sum(lambda[21:50])) else 0
n_RI[, 7, 40] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:39]), lambda[40] / sum(lambda[21:50])) else 0
n_RI[, 7, 41] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:40]), lambda[41] / sum(lambda[21:50])) else 0
n_RI[, 7, 42] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:41]), lambda[42] / sum(lambda[21:50])) else 0
n_RI[, 7, 43] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:42]), lambda[43] / sum(lambda[21:50])) else 0
n_RI[, 7, 44] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:43]), lambda[44] / sum(lambda[21:50])) else 0
n_RI[, 7, 45] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:44]), lambda[45] / sum(lambda[21:50])) else 0
n_RI[, 7, 46] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:45]), lambda[46] / sum(lambda[21:50])) else 0
n_RI[, 7, 47] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:46]), lambda[47] / sum(lambda[21:50])) else 0
n_RI[, 7, 48] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:47]), lambda[48] / sum(lambda[21:50])) else 0
n_RI[, 7, 49] <- if (sum(lambda[21:50]) > 0) Binomial(R_out[i, 7] - sum(n_RI[i, 7, 21:48]), lambda[49] / sum(lambda[21:50])) else 0
n_RI[, 7, 50] <- R_out[i, 7] - sum(n_RI[i, 7, 21:49])
n_RI[, 8, 1] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8], lambda[1] / sum(lambda[1:60])) else 0
n_RI[, 8, 2] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:1]), lambda[2] / sum(lambda[1:60])) else 0
n_RI[, 8, 3] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:2]), lambda[3] / sum(lambda[1:60])) else 0
n_RI[, 8, 4] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:3]), lambda[4] / sum(lambda[1:60])) else 0
n_RI[, 8, 5] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:4]), lambda[5] / sum(lambda[1:60])) else 0
n_RI[, 8, 6] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:5]), lambda[6] / sum(lambda[1:60])) else 0
n_RI[, 8, 7] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:6]), lambda[7] / sum(lambda[1:60])) else 0
n_RI[, 8, 8] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:7]), lambda[8] / sum(lambda[1:60])) else 0
n_RI[, 8, 9] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:8]), lambda[9] / sum(lambda[1:60])) else 0
n_RI[, 8, 10] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:9]), lambda[10] / sum(lambda[1:60])) else 0
n_RI[, 8, 11] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:10]), lambda[11] / sum(lambda[1:60])) else 0
n_RI[, 8, 12] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:11]), lambda[12] / sum(lambda[1:60])) else 0
n_RI[, 8, 13] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:12]), lambda[13] / sum(lambda[1:60])) else 0
n_RI[, 8, 14] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:13]), lambda[14] / sum(lambda[1:60])) else 0
n_RI[, 8, 15] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:14]), lambda[15] / sum(lambda[1:60])) else 0
n_RI[, 8, 16] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:15]), lambda[16] / sum(lambda[1:60])) else 0
n_RI[, 8, 17] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:16]), lambda[17] / sum(lambda[1:60])) else 0
n_RI[, 8, 18] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:17]), lambda[18] / sum(lambda[1:60])) else 0
n_RI[, 8, 19] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:18]), lambda[19] / sum(lambda[1:60])) else 0
n_RI[, 8, 20] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:19]), lambda[20] / sum(lambda[1:60])) else 0
n_RI[, 8, 51] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:20]), lambda[51] / sum(lambda[1:60])) else 0
n_RI[, 8, 52] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:51]), lambda[52] / sum(lambda[1:60])) else 0
n_RI[, 8, 53] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:52]), lambda[53] / sum(lambda[1:60])) else 0
n_RI[, 8, 54] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:53]), lambda[54] / sum(lambda[1:60])) else 0
n_RI[, 8, 55] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:54]), lambda[55] / sum(lambda[1:60])) else 0
n_RI[, 8, 56] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:55]), lambda[56] / sum(lambda[1:60])) else 0
n_RI[, 8, 57] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:56]), lambda[57] / sum(lambda[1:60])) else 0
n_RI[, 8, 58] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:57]), lambda[58] / sum(lambda[1:60])) else 0
n_RI[, 8, 59] <- if (sum(lambda[1:60]) > 0) Binomial(R_out[i, 8] - sum(n_RI[i, 8, 1:58]), lambda[59] / sum(lambda[1:60])) else 0
n_RI[, 8, 60] <- R_out[i, 8] - sum(n_RI[i, 8, 1:59])
n_RI[, 9, 1] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9], lambda[1] / sum(lambda[1:50])) else 0
n_RI[, 9, 2] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:1]), lambda[2] / sum(lambda[1:50])) else 0
n_RI[, 9, 3] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:2]), lambda[3] / sum(lambda[1:50])) else 0
n_RI[, 9, 4] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:3]), lambda[4] / sum(lambda[1:50])) else 0
n_RI[, 9, 5] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:4]), lambda[5] / sum(lambda[1:50])) else 0
n_RI[, 9, 6] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:5]), lambda[6] / sum(lambda[1:50])) else 0
n_RI[, 9, 7] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:6]), lambda[7] / sum(lambda[1:50])) else 0
n_RI[, 9, 8] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:7]), lambda[8] / sum(lambda[1:50])) else 0
n_RI[, 9, 9] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:8]), lambda[9] / sum(lambda[1:50])) else 0
n_RI[, 9, 10] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:9]), lambda[10] / sum(lambda[1:50])) else 0
n_RI[, 9, 11] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:10]), lambda[11] / sum(lambda[1:50])) else 0
n_RI[, 9, 12] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:11]), lambda[12] / sum(lambda[1:50])) else 0
n_RI[, 9, 13] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:12]), lambda[13] / sum(lambda[1:50])) else 0
n_RI[, 9, 14] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:13]), lambda[14] / sum(lambda[1:50])) else 0
n_RI[, 9, 15] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:14]), lambda[15] / sum(lambda[1:50])) else 0
n_RI[, 9, 16] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:15]), lambda[16] / sum(lambda[1:50])) else 0
n_RI[, 9, 17] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:16]), lambda[17] / sum(lambda[1:50])) else 0
n_RI[, 9, 18] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:17]), lambda[18] / sum(lambda[1:50])) else 0
n_RI[, 9, 19] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:18]), lambda[19] / sum(lambda[1:50])) else 0
n_RI[, 9, 20] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:19]), lambda[20] / sum(lambda[1:50])) else 0
n_RI[, 9, 41] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:20]), lambda[41] / sum(lambda[1:50])) else 0
n_RI[, 9, 42] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:41]), lambda[42] / sum(lambda[1:50])) else 0
n_RI[, 9, 43] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:42]), lambda[43] / sum(lambda[1:50])) else 0
n_RI[, 9, 44] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:43]), lambda[44] / sum(lambda[1:50])) else 0
n_RI[, 9, 45] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:44]), lambda[45] / sum(lambda[1:50])) else 0
n_RI[, 9, 46] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:45]), lambda[46] / sum(lambda[1:50])) else 0
n_RI[, 9, 47] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:46]), lambda[47] / sum(lambda[1:50])) else 0
n_RI[, 9, 48] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:47]), lambda[48] / sum(lambda[1:50])) else 0
n_RI[, 9, 49] <- if (sum(lambda[1:50]) > 0) Binomial(R_out[i, 9] - sum(n_RI[i, 9, 1:48]), lambda[49] / sum(lambda[1:50])) else 0
n_RI[, 9, 50] <- R_out[i, 9] - sum(n_RI[i, 9, 1:49])
n_RI[, 10, 1] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10], lambda[1] / sum(lambda[1:40])) else 0
n_RI[, 10, 2] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:1]), lambda[2] / sum(lambda[1:40])) else 0
n_RI[, 10, 3] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:2]), lambda[3] / sum(lambda[1:40])) else 0
n_RI[, 10, 4] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:3]), lambda[4] / sum(lambda[1:40])) else 0
n_RI[, 10, 5] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:4]), lambda[5] / sum(lambda[1:40])) else 0
n_RI[, 10, 6] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:5]), lambda[6] / sum(lambda[1:40])) else 0
n_RI[, 10, 7] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:6]), lambda[7] / sum(lambda[1:40])) else 0
n_RI[, 10, 8] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:7]), lambda[8] / sum(lambda[1:40])) else 0
n_RI[, 10, 9] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:8]), lambda[9] / sum(lambda[1:40])) else 0
n_RI[, 10, 10] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:9]), lambda[10] / sum(lambda[1:40])) else 0
n_RI[, 10, 11] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:10]), lambda[11] / sum(lambda[1:40])) else 0
n_RI[, 10, 12] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:11]), lambda[12] / sum(lambda[1:40])) else 0
n_RI[, 10, 13] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:12]), lambda[13] / sum(lambda[1:40])) else 0
n_RI[, 10, 14] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:13]), lambda[14] / sum(lambda[1:40])) else 0
n_RI[, 10, 15] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:14]), lambda[15] / sum(lambda[1:40])) else 0
n_RI[, 10, 16] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:15]), lambda[16] / sum(lambda[1:40])) else 0
n_RI[, 10, 17] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:16]), lambda[17] / sum(lambda[1:40])) else 0
n_RI[, 10, 18] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:17]), lambda[18] / sum(lambda[1:40])) else 0
n_RI[, 10, 19] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:18]), lambda[19] / sum(lambda[1:40])) else 0
n_RI[, 10, 20] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:19]), lambda[20] / sum(lambda[1:40])) else 0
n_RI[, 10, 21] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:20]), lambda[21] / sum(lambda[1:40])) else 0
n_RI[, 10, 22] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:21]), lambda[22] / sum(lambda[1:40])) else 0
n_RI[, 10, 23] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:22]), lambda[23] / sum(lambda[1:40])) else 0
n_RI[, 10, 24] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:23]), lambda[24] / sum(lambda[1:40])) else 0
n_RI[, 10, 25] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:24]), lambda[25] / sum(lambda[1:40])) else 0
n_RI[, 10, 26] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:25]), lambda[26] / sum(lambda[1:40])) else 0
n_RI[, 10, 27] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:26]), lambda[27] / sum(lambda[1:40])) else 0
n_RI[, 10, 28] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:27]), lambda[28] / sum(lambda[1:40])) else 0
n_RI[, 10, 29] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:28]), lambda[29] / sum(lambda[1:40])) else 0
n_RI[, 10, 30] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:29]), lambda[30] / sum(lambda[1:40])) else 0
n_RI[, 10, 31] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:30]), lambda[31] / sum(lambda[1:40])) else 0
n_RI[, 10, 32] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:31]), lambda[32] / sum(lambda[1:40])) else 0
n_RI[, 10, 33] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:32]), lambda[33] / sum(lambda[1:40])) else 0
n_RI[, 10, 34] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:33]), lambda[34] / sum(lambda[1:40])) else 0
n_RI[, 10, 35] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:34]), lambda[35] / sum(lambda[1:40])) else 0
n_RI[, 10, 36] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:35]), lambda[36] / sum(lambda[1:40])) else 0
n_RI[, 10, 37] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:36]), lambda[37] / sum(lambda[1:40])) else 0
n_RI[, 10, 38] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:37]), lambda[38] / sum(lambda[1:40])) else 0
n_RI[, 10, 39] <- if (sum(lambda[1:40]) > 0) Binomial(R_out[i, 10] - sum(n_RI[i, 10, 1:38]), lambda[39] / sum(lambda[1:40])) else 0
n_RI[, 10, 40] <- R_out[i, 10] - sum(n_RI[i, 10, 1:39])
n_RI[, 11, 51] <- if (sum(lambda[51:60]) > 0) Binomial(R_out[i, 11], lambda[51] / sum(lambda[51:60])) else 0
n_RI[, 11, 52] <- if (sum(lambda[51:60]) > 0) Binomial(R_out[i, 11] - sum(n_RI[i, 11, 51:51]), lambda[52] / sum(lambda[51:60])) else 0
n_RI[, 11, 53] <- if (sum(lambda[51:60]) > 0) Binomial(R_out[i, 11] - sum(n_RI[i, 11, 51:52]), lambda[53] / sum(lambda[51:60])) else 0
n_RI[, 11, 54] <- if (sum(lambda[51:60]) > 0) Binomial(R_out[i, 11] - sum(n_RI[i, 11, 51:53]), lambda[54] / sum(lambda[51:60])) else 0
n_RI[, 11, 55] <- if (sum(lambda[51:60]) > 0) Binomial(R_out[i, 11] - sum(n_RI[i, 11, 51:54]), lambda[55] / sum(lambda[51:60])) else 0
n_RI[, 11, 56] <- if (sum(lambda[51:60]) > 0) Binomial(R_out[i, 11] - sum(n_RI[i, 11, 51:55]), lambda[56] / sum(lambda[51:60])) else 0
n_RI[, 11, 57] <- if (sum(lambda[51:60]) > 0) Binomial(R_out[i, 11] - sum(n_RI[i, 11, 51:56]), lambda[57] / sum(lambda[51:60])) else 0
n_RI[, 11, 58] <- if (sum(lambda[51:60]) > 0) Binomial(R_out[i, 11] - sum(n_RI[i, 11, 51:57]), lambda[58] / sum(lambda[51:60])) else 0
n_RI[, 11, 59] <- if (sum(lambda[51:60]) > 0) Binomial(R_out[i, 11] - sum(n_RI[i, 11, 51:58]), lambda[59] / sum(lambda[51:60])) else 0
n_RI[, 11, 60] <- R_out[i, 11] - sum(n_RI[i, 11, 51:59])
n_RI[, 12, 41] <- if (sum(lambda[41:50]) > 0) Binomial(R_out[i, 12], lambda[41] / sum(lambda[41:50])) else 0
n_RI[, 12, 42] <- if (sum(lambda[41:50]) > 0) Binomial(R_out[i, 12] - sum(n_RI[i, 12, 41:41]), lambda[42] / sum(lambda[41:50])) else 0
n_RI[, 12, 43] <- if (sum(lambda[41:50]) > 0) Binomial(R_out[i, 12] - sum(n_RI[i, 12, 41:42]), lambda[43] / sum(lambda[41:50])) else 0
n_RI[, 12, 44] <- if (sum(lambda[41:50]) > 0) Binomial(R_out[i, 12] - sum(n_RI[i, 12, 41:43]), lambda[44] / sum(lambda[41:50])) else 0
n_RI[, 12, 45] <- if (sum(lambda[41:50]) > 0) Binomial(R_out[i, 12] - sum(n_RI[i, 12, 41:44]), lambda[45] / sum(lambda[41:50])) else 0
n_RI[, 12, 46] <- if (sum(lambda[41:50]) > 0) Binomial(R_out[i, 12] - sum(n_RI[i, 12, 41:45]), lambda[46] / sum(lambda[41:50])) else 0
n_RI[, 12, 47] <- if (sum(lambda[41:50]) > 0) Binomial(R_out[i, 12] - sum(n_RI[i, 12, 41:46]), lambda[47] / sum(lambda[41:50])) else 0
n_RI[, 12, 48] <- if (sum(lambda[41:50]) > 0) Binomial(R_out[i, 12] - sum(n_RI[i, 12, 41:47]), lambda[48] / sum(lambda[41:50])) else 0
n_RI[, 12, 49] <- if (sum(lambda[41:50]) > 0) Binomial(R_out[i, 12] - sum(n_RI[i, 12, 41:48]), lambda[49] / sum(lambda[41:50])) else 0
n_RI[, 12, 50] <- R_out[i, 12] - sum(n_RI[i, 12, 41:49])
n_RI[, 13, 21] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13], lambda[21] / sum(lambda[21:40])) else 0
n_RI[, 13, 22] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:21]), lambda[22] / sum(lambda[21:40])) else 0
n_RI[, 13, 23] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:22]), lambda[23] / sum(lambda[21:40])) else 0
n_RI[, 13, 24] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:23]), lambda[24] / sum(lambda[21:40])) else 0
n_RI[, 13, 25] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:24]), lambda[25] / sum(lambda[21:40])) else 0
n_RI[, 13, 26] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:25]), lambda[26] / sum(lambda[21:40])) else 0
n_RI[, 13, 27] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:26]), lambda[27] / sum(lambda[21:40])) else 0
n_RI[, 13, 28] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:27]), lambda[28] / sum(lambda[21:40])) else 0
n_RI[, 13, 29] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:28]), lambda[29] / sum(lambda[21:40])) else 0
n_RI[, 13, 30] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:29]), lambda[30] / sum(lambda[21:40])) else 0
n_RI[, 13, 31] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:30]), lambda[31] / sum(lambda[21:40])) else 0
n_RI[, 13, 32] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:31]), lambda[32] / sum(lambda[21:40])) else 0
n_RI[, 13, 33] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:32]), lambda[33] / sum(lambda[21:40])) else 0
n_RI[, 13, 34] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:33]), lambda[34] / sum(lambda[21:40])) else 0
n_RI[, 13, 35] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:34]), lambda[35] / sum(lambda[21:40])) else 0
n_RI[, 13, 36] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:35]), lambda[36] / sum(lambda[21:40])) else 0
n_RI[, 13, 37] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:36]), lambda[37] / sum(lambda[21:40])) else 0
n_RI[, 13, 38] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:37]), lambda[38] / sum(lambda[21:40])) else 0
n_RI[, 13, 39] <- if (sum(lambda[21:40]) > 0) Binomial(R_out[i, 13] - sum(n_RI[i, 13, 21:38]), lambda[39] / sum(lambda[21:40])) else 0
n_RI[, 13, 40] <- R_out[i, 13] - sum(n_RI[i, 13, 21:39])
n_RI[, 14, 1] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14], lambda[1] / sum(lambda[1:20])) else 0
n_RI[, 14, 2] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:1]), lambda[2] / sum(lambda[1:20])) else 0
n_RI[, 14, 3] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:2]), lambda[3] / sum(lambda[1:20])) else 0
n_RI[, 14, 4] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:3]), lambda[4] / sum(lambda[1:20])) else 0
n_RI[, 14, 5] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:4]), lambda[5] / sum(lambda[1:20])) else 0
n_RI[, 14, 6] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:5]), lambda[6] / sum(lambda[1:20])) else 0
n_RI[, 14, 7] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:6]), lambda[7] / sum(lambda[1:20])) else 0
n_RI[, 14, 8] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:7]), lambda[8] / sum(lambda[1:20])) else 0
n_RI[, 14, 9] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:8]), lambda[9] / sum(lambda[1:20])) else 0
n_RI[, 14, 10] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:9]), lambda[10] / sum(lambda[1:20])) else 0
n_RI[, 14, 11] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:10]), lambda[11] / sum(lambda[1:20])) else 0
n_RI[, 14, 12] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:11]), lambda[12] / sum(lambda[1:20])) else 0
n_RI[, 14, 13] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:12]), lambda[13] / sum(lambda[1:20])) else 0
n_RI[, 14, 14] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:13]), lambda[14] / sum(lambda[1:20])) else 0
n_RI[, 14, 15] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:14]), lambda[15] / sum(lambda[1:20])) else 0
n_RI[, 14, 16] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:15]), lambda[16] / sum(lambda[1:20])) else 0
n_RI[, 14, 17] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:16]), lambda[17] / sum(lambda[1:20])) else 0
n_RI[, 14, 18] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:17]), lambda[18] / sum(lambda[1:20])) else 0
n_RI[, 14, 19] <- if (sum(lambda[1:20]) > 0) Binomial(R_out[i, 14] - sum(n_RI[i, 14, 1:18]), lambda[19] / sum(lambda[1:20])) else 0
n_RI[, 14, 20] <- R_out[i, 14] - sum(n_RI[i, 14, 1:19])
