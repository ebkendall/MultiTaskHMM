# ------------------------------------------------------------------------------
# Simulation setup: ------------------------------------------------------------
# Assume that we have 3 states, the response is of dimension 5, and the 
# covariance matrices for each state are diagonal (for now)
# ------------------------------------------------------------------------------

library(mvtnorm)
library(LaplacesDemon)
set.seed(2023)

# Number of datasets to simulate
n_sim = 100

# Number of environments to simulate
n_env = 3

# Sample size
N = 100

# Number of states
S = 3

# *************** Environment 1 ***************
# True parameter values
true_par_e1 = c(c(matrix(c(-2, -1.5,
                           -2,  1.0,
                           -2,  1.5,
                           -2,  1.5,
                           -2,  0.5,
                           -2, -2.0), ncol = 2, byrow = T)),
             0.45, 0.1, 0.45,                   # initial probabilities (same)
             0, 0, 0, 0, 0,                     # mu_1 (different)
             1, 1, 1, 1, 1,                     # mu_2 (different)
             2, 2, 2, 2, 2,                     # mu_3 (different)
             c(diag(5)),                        # Sigma_1 (different)
             c(diag(5)),                        # Sigma_2 (different)
             c(diag(5)))                        # Sigma_3 (different)

# *************** Environment 2 ***************
# True parameter values
true_par_e2 = c(c(matrix(c(-2, -1.5,
                           -2,  1.0,
                           -2,  1.5,
                           -2,  1.5,
                           -2,  0.5,
                           -2, -2.0), ncol = 2, byrow = T)),
                0.4, 0.3, 0.3,                  # initial probabilities (same)
                -1, -1, -1, -1, -1,             # mu_1 (different)  
                1, 1, 1, 1, 1,                  # mu_2 (different)
                -2, -2, -2, -2, -2,             # mu_3 (different)
                c(diag(rep(2,5))),              # Sigma_1 (different)
                c(diag(rep(0.5, 5))),           # Sigma_2 (different)
                c(diag(rep(1.5, 5))))           # Sigma_3 (different)

# *************** Environment 3 ***************
# True parameter values
true_par_e3 = c(c(matrix(c(-2, -1.5,
                           -2,  1.0,
                           -2,  1.5,
                           -2,  1.5,
                           -2,  0.5,
                           -2, -2.0), ncol = 2, byrow = T)),
                0.4, 0.35, 0.25,                # initial probabilities (same)
                0.5, 0.5, 0.5, 0.5, 0.5,        # mu_1 (different) 
                -1, -1, -1, -1, -1,             # mu_2 (different)
                2, 2, 2, 2, 2,                  # mu_3 (different)
                c(rinvwishart(9, diag(5))),     # Sigma_1 (different)
                c(rinvwishart(9, diag(5))),     # Sigma_2 (different)
                c(rinvwishart(9, diag(5))))     # Sigma_3 (different)

# *************** All environments ***************
true_par = list(true_par_e1, true_par_e2, true_par_e3)
save(true_par, file = 'Data/true_par.rda')


# Indexing of the parameter vector
par_index = list(t_p = 1:12, init_pi = 13:15, 
                 mu_1 = 16:20, mu_2 = 21:25, mu_3 = 26:30, 
                 Sig_1 = 31:55, Sig_2 = 56:80, Sig_3 = 81:105)

P_mat <- function(zeta, x) {
    zetaMat = matrix(zeta, ncol = 2)
    x_mat = matrix(c(1, x), nrow = 2)
    
    p_comp = c(zetaMat %*% x_mat)
    p_comp = exp(p_comp)
    
    P = matrix(c(        1, p_comp[1], p_comp[2],
                 p_comp[3],         1, p_comp[4],
                 p_comp[5], p_comp[6],         1), ncol = 3, byrow = T)
    
    P = P / rowSums(P)
    
    return(P)
}

for(seed in 1:n_sim) {

    set.seed(seed)

    for(e in 1:n_env) {
        
        sim_data = NULL
        
        indiv_sample_size = rpois(n = N, lambda = 15)
        
        Sigma_1 = matrix(true_par[[e]][par_index$Sig_1], nrow = 5)
        Sigma_2 = matrix(true_par[[e]][par_index$Sig_2], nrow = 5)
        Sigma_3 = matrix(true_par[[e]][par_index$Sig_3], nrow = 5)
        
        # Initial State probabilities
        init_prob = true_par[[e]][par_index$init_pi]
        
        for(i in 1:N) {
            
            # Transition probability matrix
            zeta = true_par[[1]][par_index$t_p]
            
            # Covariate distribution will vary between environments
            if(e == 1) {
                x_i = runif(1, min = 1, max = 2)
            } else if(e == 2) {
                x_i = runif(1, min = -2, max = -1)
            } else {
                x_i = runif(1, min = -2, max = 2)
            }
            
            P = P_mat(zeta, x_i)
            
            n_i = indiv_sample_size[i]
            
            # Simulate the latent state sequence ----------------------------------
            # First the initial state
            s_i = sample(1:3, size=1, prob=init_prob)
            
            # Second the remaining states
            for(j in 2:n_i) {
                s_i = c( s_i, sample(1:3, size=1, prob=P[tail(s_i,1),]))
            }
            # ---------------------------------------------------------------------
            
            # Simulate the response sequence --------------------------------------
            Y_i = matrix(nrow = n_i, ncol = 5)
            
            for(j in 1:n_i) {
                if(s_i[j] == 1) {
                    Y_i[j, ] = rmvnorm(n=1, mean = true_par[[e]][par_index$mu_1], 
                                       sigma = Sigma_1)
                } else if(s_i[j] == 2) {
                    Y_i[j, ] = rmvnorm(n=1, mean = true_par[[e]][par_index$mu_2], 
                                       sigma = Sigma_2)
                } else {
                    Y_i[j, ] = rmvnorm(n=1, mean = true_par[[e]][par_index$mu_3], 
                                       sigma = Sigma_3)
                }
            }
            # ---------------------------------------------------------------------
            
            sim_data_sub = cbind(rep(i, n_i), Y_i, s_i, rep(x_i, n_i))
            
            sim_data = rbind(sim_data, sim_data_sub)
            
        }
        
        colnames(sim_data) = c("id", "y1", "y2", "y3", "y4", "y5", "state", "x")
        
        save(sim_data, file = paste0("Data/sim_data_", seed, "_e", e, ".rda"))
        
        print(paste0("Environment: ", e, " -------------------------"))
        cat('\n','Proption of occurances in each state:','\n')
        print(table(sim_data[,'state']))
        print(table(sim_data[,'state'])/dim(sim_data)[1])
        cat('\n')
        
        count_transitions = matrix(0, nrow=3, ncol=3)
        total_trans = 0
        for(i in unique(sim_data[,"id"])){
            
            b_i_mle = as.numeric(c(sim_data[sim_data[,"id"] == i, 'state']))
            
            for(t in 1:(length(b_i_mle) - 1)) {
                count_transitions[b_i_mle[t], b_i_mle[t+1]] = 
                    count_transitions[b_i_mle[t], b_i_mle[t+1]] + 1
                total_trans = total_trans + 1
            }
        }
        
        print("Transition distribution")
        print(count_transitions)
    }
}
