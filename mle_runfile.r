source("mle_routine.r")

args = commandArgs(TRUE)
ind = as.numeric(args[1])

set.seed(ind)
print(ind)

simulation = T

# Initial values
init_val = c(1/4, 1/4, 1/4, 1/4, 1/4, 1/4, # transition probs
             1/6, 1/6,                     # initial probabilities
             0, 0, 0, 0, 0,                # mu_1
             1, 1, 1, 1, 1,                # mu_2
             2, 2, 2, 2, 2,                # mu_3
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

start_t_big = Sys.time()
par_est <- baum_welch_one_environment(init_val, par_index, y, id)
end_t_big = Sys.time(); print(paste0("Total time: ", end_t_big - start_t_big))

print("Estimated P")
prob_val = par_est[par_index$t_p]
P = matrix(c(1 - prob_val[1] - prob_val[2], prob_val[1], prob_val[2],
             prob_val[3], 1 - prob_val[3] - prob_val[4], prob_val[4],
             prob_val[5], prob_val[6], 1 - prob_val[5] - prob_val[6]),
           nrow = 3, byrow = T)
print(P)

print("Estimated initial prob.")
init_val = par_est[par_index$init_pi]
init = c(1 - sum(init_val), init_val[1], init_val[2])
print(init)

print("Estimated mu_1")
print(par_est[par_index$mu_1])
print("Estimated mu_2")
print(par_est[par_index$mu_2])
print("Estimated mu_3")
print(par_est[par_index$mu_3])

print("Estimated Sigma 1")
print(matrix(par_est[par_index$Sig_1], nrow = 5))
print("Estimated Sigma 2")
print(matrix(par_est[par_index$Sig_2], nrow = 5))
print("Estimated Sigma 3")
print(matrix(par_est[par_index$Sig_3], nrow = 5))