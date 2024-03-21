source("mle_routine.r")

args = commandArgs(TRUE)
ind = as.numeric(args[1])

set.seed(ind)
print(ind)

simulation = T

# Initial values
init_val = c(-2, -2, -2, -2, -2, -2, # logit transition probs
             -2, -2,                 # logit initial probabilities
              0, 0, 0, 0, 0,         # mu_1
              2, 2, 2, 2, 2,         # mu_2
              4, 4, 4, 4, 4,         # mu_3
              0, 0, 0, 0, 0,         # log_sigma_1
              0, 0, 0, 0, 0,         # log_sigma_2
              0, 0, 0, 0, 0)         # log_sigma_3

# Indexing of the parameter vector
par_index = list(logit_t_p = 1:6, logit_pi = 7:8, 
                 mu_1 = 9:13, mu_2 = 14:18, mu_3 = 19:23, 
                 log_sig_1 = 24:28, log_sig_2 = 29:33, log_sig_3 = 34:38)

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