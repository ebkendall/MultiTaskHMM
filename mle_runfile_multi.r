source("mle_routine_multi.r")

args = commandArgs(TRUE)
ind = as.numeric(args[1])

n_env = 1

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

init_val_list <- list()
for(e in 1:n_env) init_val_list[[e]] = init_val

# Indexing of the parameter vector
par_index = list(t_p = 1:6, init_pi = 7:8, 
                 mu_1 = 9:13, mu_2 = 14:18, mu_3 = 19:23, 
                 Sig_1 = 24:48, Sig_2 = 49:73, Sig_3 = 74:98)

id = y = list()

if(simulation) {
    for(e in 1:n_env) {
        # Load simulated data
        load(paste0('Data/sim_data_', ind, "_e", e, '.rda'))
        
        id[[e]] = sim_data[,'id']
        y[[e]]  = sim_data[,c('y1', 'y2', 'y3', 'y4', 'y5')]   
    }
} else {
    # Load whatever real data we get
}

start_t_big = Sys.time()
par = init_val_list
par_est <- baum_welch_multi_environment(init_val_list, par_index, y, id, n_env)
end_t_big = Sys.time(); print(paste0("Total time: ", end_t_big - start_t_big))

model_out = list(par_est, end_t_big - start_t_big, par_index)

save(model_out, file = paste0('Model_out/model_out_multi_e', n_env, "_ind", ind, '.rda'))

for(e in 1:n_env) {
    print(paste0("Environment ", e))
    par_est_e = par_est[[e]]

    print("Estimated P")
    prob_val = par_est_e[par_index$t_p]
    P = matrix(c(1 - prob_val[1] - prob_val[2], prob_val[1], prob_val[2],
                prob_val[3], 1 - prob_val[3] - prob_val[4], prob_val[4],
                prob_val[5], prob_val[6], 1 - prob_val[5] - prob_val[6]),
            nrow = 3, byrow = T)
    print(P)

    print("Estimated initial prob.")
    init_val = par_est_e[par_index$init_pi]
    init = c(1 - sum(init_val), init_val[1], init_val[2])
    print(init)

    print("Estimated mu_1")
    print(par_est_e[par_index$mu_1])
    print("Estimated mu_2")
    print(par_est_e[par_index$mu_2])
    print("Estimated mu_3")
    print(par_est_e[par_index$mu_3])

    print("Estimated Sigma 1")
    print(matrix(par_est_e[par_index$Sig_1], nrow = 5))
    print("Estimated Sigma 2")
    print(matrix(par_est_e[par_index$Sig_2], nrow = 5))
    print("Estimated Sigma 3")
    print(matrix(par_est_e[par_index$Sig_3], nrow = 5))
}