source("mle_routine.r")

args = commandArgs(TRUE)
ind = as.numeric(args[1])

set.seed(ind)
print(ind)

simulation = T

# Initial values
init_val = c(1/3, 1/3, 1/3, 1/3, 1/3, 1/3, # transition probs
             1/3, 1/3,                     # initial probabilities
             0, 0, 0, 0, 0,                # mu_1
             2, 2, 2, 2, 2,                # mu_2
             4, 4, 4, 4, 4,                # mu_3
             c(diag(5)),                   # Sigma_1
             c(diag(5)),                   # Sigma_2
             c(diag(5)))                   # Sigma_3

# Indexing of the parameter vector
par_index = list(t_p = 1:6, init_pi = 7:8, 
                 mu_1 = 9:13, mu_2 = 14:18, mu_3 = 19:23, 
                 Sig_1 = 24:48, Sig_2 = 49:73, Sig_3 = 74:98)

id = y = NULL

if(simulation) {
    # Load simulated data
    load(paste0('Data/sim_data_', ind, '.rda'))

    id = sim_data[,'id']
    y  = sim_data[,c('y1', 'y2', 'y3', 'y4', 'y5')]
} else {
    # Load whatever real data we get
}

par_est <- baum_welch_one_environment(init_val, par_index, y, id)