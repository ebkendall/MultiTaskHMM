# ------------------------------------------------------------------------------
# Simulation setup: ------------------------------------------------------------
# Assume that we have 3 states, the response is of dimension 5, and the 
# covariance matrices for each state are diagonal (for now)
# ------------------------------------------------------------------------------

library(mvtnorm)

# Number of datasets to simulate
n_sim = 1

# Sample size
N = 100

# Number of states
S = 3

# True parameter values
true_par = c(-2, -2, -2, -2, -2, -2, # logit transition probs
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

# Transition probability matrix
# *** NOTE: this will be more interesting is this is in terms of covariates
q1   = exp(true_par[par_index$logit_t_p][1])
q2   = exp(true_par[par_index$logit_t_p][2])
q3   = exp(true_par[par_index$logit_t_p][3])
q4   = exp(true_par[par_index$logit_t_p][4])
q5   = exp(true_par[par_index$logit_t_p][5])
q6   = exp(true_par[par_index$logit_t_p][6])

Q = matrix(c(1,  q1, q2,
            q3,   1, q4,
            q5,  q6,  1), ncol=3, byrow=T)
P = Q / rowSums(Q)

# Initial State probabilities
init_prob = c(1, exp(true_par[par_index$logit_pi][1]), 
                 exp(true_par[par_index$logit_pi][2]))
init_prob = init_prob / sum(init_prob)

# Response Covariance matrices
Sigma_1 = diag(exp(true_par[par_index$log_sig_1]))
Sigma_2 = diag(exp(true_par[par_index$log_sig_2]))
Sigma_3 = diag(exp(true_par[par_index$log_sig_3]))

for(seed in 1:n_sim) {

    set.seed(seed)

    sim_data = NULL

    indiv_sample_size = rpois(n = N, lambda = 15)

    for(i in 1:N) {

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
                Y_i[j, ] = rmvnorm(n=1, mean = true_par[par_index$mu_1], 
                                        sigma = Sigma_1)
            } else if(s_i[j] == 2) {
                Y_i[j, ] = rmvnorm(n=1, mean = true_par[par_index$mu_2], 
                                        sigma = Sigma_2)
            } else {
                Y_i[j, ] = rmvnorm(n=1, mean = true_par[par_index$mu_3], 
                                        sigma = Sigma_3)
            }
        }
        # ---------------------------------------------------------------------

        sim_data_sub = cbind(rep(i, n_i), Y_i, s_i)

        sim_data = rbind(sim_data, sim_data_sub)

    }

    colnames(sim_data) = c("id", "y1", "y2", "y3", "y4", "y5", "state")

    save(sim_data, file = paste0("Data/sim_data_", seed, ".rda"))

    cat('\n','Proption of occurances in each state:','\n')
    print(table(sim_data[,'state']))
    print(table(sim_data[,'state'])/dim(sim_data)[1])
    cat('\n')
}
