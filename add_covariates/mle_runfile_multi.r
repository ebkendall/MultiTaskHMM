source("mle_routine_multi.r")

# args = commandArgs(TRUE)
# ind = as.numeric(args[1])

ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

n_env = 3

set.seed(ind)
print(ind)

simulation = T

# Initial values
init_val = c(c(matrix(c(-2, -1.5,
                        -2,  1.0,
                        -2,  1.5,
                        -2,  1.5,
                        -2,  0.5,
                        -2, -2.0), ncol = 2, byrow = T)),
             0.40, 0.30, 0.30,             # initial probabilities
             0, 0, 0, 0, 0,                # mu_1
             1, 1, 1, 1, 1,                # mu_2
             2, 2, 2, 2, 2,                # mu_3
             c(diag(5)),                   # Sigma_1
             c(diag(5)),                   # Sigma_2
             c(diag(5)))                   # Sigma_3

init_val_list <- list()
for(e in 1:n_env) init_val_list[[e]] = init_val

# Indexing of the parameter vector
par_index = list(t_p = 1:12, init_pi = 13:15, 
                 mu_1 = 16:20, mu_2 = 21:25, mu_3 = 26:30, 
                 Sig_1 = 31:55, Sig_2 = 56:80, Sig_3 = 81:105)

id = y = x = list()

if(simulation) {
    for(e in 1:n_env) {
        # Load simulated data
        load(paste0('Data/sim_data_', ind, "_e", e, '.rda'))
        
        id[[e]] = sim_data[,'id']
        y[[e]]  = sim_data[,c('y1', 'y2', 'y3', 'y4', 'y5')]  
        x[[e]]  = cbind(1, sim_data[,'x'])
    }
} else {
    # Load whatever real data we get
}

start_t_big = Sys.time()
par = init_val_list
par_est <- baum_welch_multi_environment(init_val_list, par_index, y, id, x, n_env)
end_t_big = Sys.time(); print(paste0("Total time: ", end_t_big - start_t_big))

model_out = list(par_est, end_t_big - start_t_big, par_index)

save(model_out, file = paste0('Model_out/model_out_multi_e', n_env, "_ind", ind, '.rda'))

for(e in 1:n_env) {
    print(paste0("Environment ", e))
    par_est_e = par_est[[e]]

    print("Estimated P")
    prob_val = par_est_e[par_index$t_p]
    P = matrix(prob_val, nrow = 3)
    print(P)

    print("Estimated initial prob.")
    init_val = par_est_e[par_index$init_pi]
    init = init_val
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