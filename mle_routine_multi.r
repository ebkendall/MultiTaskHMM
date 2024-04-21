library(mvtnorm, quietly=T) 
library(foreach, quietly=T) 
library(doParallel, quietly=T)

baum_welch_multi_environment <- function(par, par_index, y, id, n_env) {
    
    eps <- 1e-4
    
    start_t = Sys.time()
    omega_k_list <- list()
    for(e in 1:n_env) {
        omega_k_list[[e]] <- omega_k_calc(par[[e]], par_index, y[[e]], id[[e]])   
    }
    end_t = Sys.time(); print(end_t - start_t)
    
    omega_k = omega_k_1 = rep(0, n_env)
    big_gamma = big_xi = list()
    for(e in 1:n_env) {
        omega_k[e] = omega_k_list[[e]][[1]]
        big_gamma[[e]] = omega_k_list[[e]][[2]]
        big_xi[[e]]    = omega_k_list[[e]][[3]]
    }
    
    loop_cont = T
    
    mpi = list(1,2,3,4,5,6,7,8, 
               c(par_index$mu_1), c(par_index$mu_2), c(par_index$mu_3),
               c(par_index$Sig_1), c(par_index$Sig_2), c(par_index$Sig_3))
    
    print("Initial omega:")
    print(omega_k)
    
    it_count = 0
    
    while(loop_cont) {
        omega_k_1 = omega_k

        for(e in 1:n_env) {
            for(j in 1:length(mpi)) {
                it_count = it_count + 1
                
                ind_j = mpi[[j]]
                
                if(sum(ind_j %in% par_index$t_p) == length(ind_j)) {
                    # transition prob. (shared across environment)
                    par[[1]][ind_j] = par[[2]][ind_j] = par[[3]][ind_j] = A_sm_update(ind_j, big_gamma, big_xi, id, n_env)
                } else if(sum(ind_j %in% par_index$init_pi) == length(ind_j)) {
                    # initial state prob. (shared across environment)
                    par[[1]][ind_j] = par[[2]][ind_j] = par[[3]][ind_j] = pi_s_update(ind_j - length(par_index$t_p) + 1, big_gamma, id, n_env)
                } else if(sum(ind_j %in% par_index$mu_1) == length(ind_j)) {
                    # mu_1
                    par[[e]][ind_j] = mu_s_update(1, big_gamma[[e]], 
                                                  y[[e]], id[[e]])
                } else if(sum(ind_j %in% par_index$mu_2) == length(ind_j)) {
                    # mu_2
                    par[[e]][ind_j] = mu_s_update(2, big_gamma[[e]], 
                                                  y[[e]], id[[e]])
                } else if(sum(ind_j %in% par_index$mu_3) == length(ind_j)) {
                    # mu_3
                    par[[e]][ind_j] = mu_s_update(3, big_gamma[[e]], 
                                                  y[[e]], id[[e]])
                } else if(sum(ind_j %in% par_index$Sig_1) == length(ind_j)) {
                    # Sig_1
                    par[[e]][ind_j] = Sigma_s_update(1, big_gamma[[e]], par[[e]], 
                                                     par_index, y[[e]], id[[e]])
                } else if(sum(ind_j %in% par_index$Sig_2) == length(ind_j)) {
                    # Sig_2
                    par[[e]][ind_j] = Sigma_s_update(2, big_gamma[[e]], par[[e]], 
                                                     par_index, y[[e]], id[[e]])
                } else {
                    # Sig_3
                    par[[e]][ind_j] = Sigma_s_update(3, big_gamma[[e]], par[[e]], 
                                                     par_index, y[[e]], id[[e]])
                }
                
                omega_k_list[[e]] <- omega_k_calc(par[[e]], par_index, y[[e]], id[[e]])
                
                omega_k[e]     = omega_k_list[[e]][[1]]
                big_gamma[[e]] = omega_k_list[[e]][[2]]
                big_xi[[e]]    = omega_k_list[[e]][[3]]
                
                # Omega needs to be updated for all env. if P or pi are updated
                if((sum(ind_j %in% par_index$t_p) == length(ind_j)) | 
                   (sum(ind_j %in% par_index$init_pi) == length(ind_j))) {
                    for(ee in 1:n_env) {
                        if(ee != e) {
                            omega_k_list[[ee]] <- omega_k_calc(par[[ee]], par_index, y[[ee]], id[[ee]]) 
                            
                            omega_k[ee] = omega_k_list[[ee]][[1]]
                            big_gamma[[ee]] = omega_k_list[[ee]][[2]]
                            big_xi[[ee]]    = omega_k_list[[ee]][[3]]
                        }
                    }
                }
                
                print(paste0("it: ", it_count, ", j: ", j))
            }   
        }

        print("prev")
        print(omega_k_1)
        print("curr")
        print(omega_k)

        if(abs(sum(omega_k) - sum(omega_k_1)) < eps) {
            loop_cont = F
            break
        }
    }
    
    return(par)
}

# Learn separately
mu_s_update <- function(s, big_gamma, y, id) {
    
    id_unique <- unique(id)
    
    mu_hat_num <- rep(0, ncol(y))
    mu_hat_den <- 0
    
    for(i in 1:length(id_unique)) {
        y_i = y[id == id_unique[i], ]
        gamma_i_s = big_gamma[[i]][,s]
        
        gamma_y = gamma_i_s * y_i
        
        mu_hat_num = mu_hat_num + colSums(gamma_y)
        mu_hat_den = mu_hat_den + sum(gamma_i_s)
    }
    
    mu_hat <- mu_hat_num / mu_hat_den
    
    return(mu_hat)
}

# Learn separately 
Sigma_s_update <- function(s, big_gamma, par, par_index, y, id) {
    
    # Initialize key components -----------------------------------------------
    id_unique <- unique(id)
    
    m_list = vector(mode = 'list', length = 3)
    m_list[[1]] = par[par_index$mu_1]
    m_list[[2]] = par[par_index$mu_2]
    m_list[[3]] = par[par_index$mu_3]
    # -------------------------------------------------------------------------
    
    sigma_hat_num <- matrix(0, nrow = 5, ncol = 5)
    sigma_hat_den <- 0
    
    for(i in 1:length(id_unique)) {
        y_i = y[id == id_unique[i], ]
        n_i = nrow(y_i)
        
        for(t in 1:n_i) {
            gamma_s_i_t <- big_gamma[[i]][t,s]
            
            y_mu_diff <- matrix(c(y_i[t, ] - m_list[[s]]), ncol = 1)
            outer_prod <- y_mu_diff %*% t(y_mu_diff)
            
            sigma_hat_num <- sigma_hat_num + gamma_s_i_t * outer_prod
            sigma_hat_den <- sigma_hat_den + gamma_s_i_t
        }
    }
    
    Sigma_hat <- sigma_hat_num / sigma_hat_den
    
    return(Sigma_hat)
}

# Pool the data
pi_s_update <- function(s, big_gamma, id, n_env) {
    
    
    pi_hat_num <- 0
    pi_hat_den <- 0
    
    for(ee in 1:n_env) {
        id_unique <- unique(id[[ee]])
        
        for(i in 1:length(id_unique)) {
            gamma_i = big_gamma[[ee]][[i]]
            
            pi_hat_num <- pi_hat_num + gamma_i[1,s]
            pi_hat_den <- pi_hat_den + sum(gamma_i[1,])
        }   
    }
    
    pi_hat <- pi_hat_num / pi_hat_den
    
    return(pi_hat)
}

# Pool the data
A_sm_update <- function(ind_j, big_gamma, big_xi, id, n_env) {
    
    pos_matrix = matrix(c(0,1,2,3,0,4,5,6,0),nrow = 3, byrow = T)
    ind_j_pos = which(pos_matrix == ind_j, arr.ind = T)
    s = ind_j_pos[1,1]
    m = ind_j_pos[1,2]
    
    a_sm_num <- 0
    a_sm_den <- 0
    
    for(ee in 1:n_env) {
        id_unique <- unique(id[[ee]])
        for(i in 1:length(id_unique)) {
            
            gamma_i = big_gamma[[ee]][[i]]
            xi_i = big_xi[[ee]][[i]]
            
            a_sm_num <- a_sm_num + sum(xi_i[[s]][[m]])
            a_sm_den <- a_sm_den + sum(gamma_i[1:(nrow(gamma_i) - 1), s])
            
        }   
    }
    
    A_sm_hat <- a_sm_num / a_sm_den
    
    return(A_sm_hat)
}

omega_k_calc <- function(par, par_index, y, id, e) {
    id_unique <- unique(id)
    
    pi_comp   <- 0
    A_comp    <- 0
    like_comp <- 0
    
    m_list = vector(mode = 'list', length = 3)
    m_list[[1]] = par[par_index$mu_1]
    m_list[[2]] = par[par_index$mu_2]
    m_list[[3]] = par[par_index$mu_3]
    
    cov_list = vector(mode = 'list', length = 3)
    cov_list[[1]] = matrix(par[par_index$Sig_1], nrow = 5)
    cov_list[[2]] = matrix(par[par_index$Sig_2], nrow = 5)
    cov_list[[3]] = matrix(par[par_index$Sig_3], nrow = 5)
    
    prob_val = par[par_index$t_p]
    P = matrix(c(1 - prob_val[1] - prob_val[2], prob_val[1], prob_val[2],
                 prob_val[3], 1 - prob_val[3] - prob_val[4], prob_val[4],
                 prob_val[5], prob_val[6], 1 - prob_val[5] - prob_val[6]),
               nrow = 3, byrow = T)
    
    init_val = par[par_index$init_pi]
    init = c(1 - sum(init_val), init_val[1], init_val[2])
    
    # After each calculation of omega, we can save these calculations for 
    # updating the other parameters -------------------------------------
    big_gamma = vector(mode = 'list', length = length(id_unique))
    big_xi    = vector(mode = 'list', length = length(id_unique))
    
    for(i in 1:length(id_unique)) {
        y_i = y[id == id_unique[i], ]
        n_i = nrow(y_i)
        
        # Calculate forward proc. & backward proc.
        alpha_mat = forward_proc_it(m_list, cov_list, init, P, y_i)
        beta_mat  = backward_proc_it(m_list, cov_list, init, P, y_i)
        
        # Calculate gamma for all time points
        gamma_mat = gamma_calc(alpha_mat, beta_mat, m_list, cov_list, init, P, y_i)
        big_gamma[[i]] = gamma_mat
        
        # pi calculation
        pi_comp = sum(gamma_mat[1,] * log(init))
        
        # transition prob calculation
        big_xi[[i]] = vector(mode = 'list', length = length(m_list))
        for(l in 1:length(m_list)) {
            big_xi[[i]][[l]] = vector(mode = 'list', length = length(m_list))
            for(j in 1:length(m_list)) {
                # The diagonal components are functions of the others
                if(j != l) {
                    xi_time = xi_calc(l, j, alpha_mat, beta_mat, m_list, cov_list, init, P, y_i)
                    big_xi[[i]][[l]][[j]] = xi_time
                    
                    A_comp = A_comp + sum(xi_time * log(P[l, j]))
                }
            }
        }
        
        # likelihood calculation
        for(t in 1:n_i) {
            for(l in 1:3) {
                like_comp = like_comp + gamma_mat[t, l] * dmvnorm(y_i[t, ], 
                                                                  mean = m_list[[l]], 
                                                                  sigma = cov_list[[l]], log = T)
            }
        }
    }
    
    Q_k = pi_comp + A_comp + like_comp
    
    omega_list = vector(mode = 'list', length = 3)
    omega_list[[1]] = Q_k
    omega_list[[2]] = big_gamma
    omega_list[[3]] = big_xi
    
    return(omega_list)
    
}

gamma_calc <- function(alpha_mat, beta_mat, m_list, cov_list, init, P, y_i) {
    
    gamma_mat = matrix(nrow = nrow(y_i), ncol = length(m_list))
    
    alpha_beta_prod = alpha_mat * beta_mat
    alpha_beta_sum = rowSums(alpha_beta_prod)
    
    for(l in 1:ncol(gamma_mat)) {
        alpha_l = alpha_mat[,l]
        beta_l  = beta_mat[,l]
        
        gamma_mat[,l] = (alpha_l * beta_l) / alpha_beta_sum
    }
    
    return(gamma_mat)
}

xi_calc <- function(l, j, alpha_mat, beta_mat, m_list, cov_list, init, P, y_i) {
    
    xi_vec = rep(0, nrow(y_i) - 1)
    
    for(t in 2:nrow(y_i)) {
        
        b_t = rep(0, length(m_list))
        for(i in 1:length(m_list)) {
            b_t[i] = dmvnorm(y_i[t, ], mean = m_list[[i]], sigma = cov_list[[i]])   
        }
        
        xi_numerator   = alpha_mat[t-1,l] * P[l,j] * beta_mat[t,j] * b_t[j]
        xi_denominator = 0
        
        for(k in 1:3) {
            for(w in 1:3) {
                xi_denominator = xi_denominator + alpha_mat[t-1, k] * P[k,w] * beta_mat[t,w] * b_t[w]
            }
        }
        
        xi_vec[t-1] = xi_numerator / xi_denominator
    }
    
    return(xi_vec)
}

forward_proc_rec <- function(t, l, m_list, cov_list, init, P, y_i) {
    
    # *** SUPER COMPUTATIONALLY EXPENSIVE *** #
    
    # m_list[[l]] = mean for state l
    # cov_list[[l]] = covariance for state l
    
    if(t == 1) {
        # Basecase
        b_l_t = dmvnorm(y_i[t, ], mean = m_list[[l]], sigma = cov_list[[l]])
        alpha_l_t = init[l] * b_l_t
        
        return(alpha_l_t)
        
    } else {
        # Recursive step
        b_l_t = dmvnorm(y_i[t, ], mean = m_list[[l]], sigma = cov_list[[l]])
        
        alpha_sum = 0
        for(j in 1:3) {
            alpha_sum = alpha_sum + P[l,j] * forward_proc_rec(t-1, j, m_list, 
                                                              cov_list, init, P,
                                                              y_i)
        }
        
        return(b_l_t * alpha_sum)
    }
}

backward_proc_rec <- function(t, l, m_list, cov_list, init, P, y_i) {
    
    # *** SUPER COMPUTATIONALLY EXPENSIVE *** #
    
    # m_list[[l]] = mean for state l
    # cov_list[[l]] = covariance for state l
    
    if(t == nrow(y_i)) {
        # Basecase
        return(1)
        
    } else {
        # Recursive step
        beta_sum = 0
        for(j in 1:3) {
            b_j_t = dmvnorm(y_i[t+1, ], mean = m_list[[j]], sigma = cov_list[[j]])
            
            beta_sum = beta_sum + P[l,j] * b_j_t * backward_proc_rec(t+1, j, 
                                                                     m_list, 
                                                                     cov_list, 
                                                                     init, P, y_i)
        }
        
        return(beta_sum)
        
    }
}

forward_proc_it <- function(m_list, cov_list, init, P, y_i) {
    
    # m_list[[l]] = mean for state l
    # cov_list[[l]] = covariance for state l
    
    alpha_mat = matrix(nrow = nrow(y_i), ncol = length(m_list))
    
    for(t in 1:nrow(y_i)) {
        if(t == 1) {
            for(l in 1:length(m_list)) {
                alpha_mat[t,l] = init[l] * dmvnorm(y_i[t, ], mean = m_list[[l]],
                                                   sigma = cov_list[[l]])
            }
        } else {
            for(l in 1:length(m_list)) {
                a_ji = c(P[,l])
                alpha_t_1 = alpha_mat[t-1, ]
                
                alpha_mat[t,l] = dmvnorm(y_i[t, ], mean = m_list[[l]], 
                                         sigma = cov_list[[l]]) * sum(a_ji * alpha_t_1)
            }
        }
    }
    return(alpha_mat)
}

backward_proc_it <- function(m_list, cov_list, init, P, y_i) {
    
    # m_list[[l]] = mean for state l
    # cov_list[[l]] = covariance for state l
    
    beta_mat = matrix(nrow = nrow(y_i), ncol = length(m_list))
    
    for(t in nrow(y_i):1) {
        if(t == nrow(y_i)) {
            beta_mat[t, ] = rep(1, ncol(beta_mat))
        } else {
            for(l in 1:length(m_list)) {
                b_j = rep(0,length(m_list))
                for(j in 1:length(m_list)) {
                    beta_t_1 = beta_mat[t+1,j]
                    a_ij = c(P[l,j])
                    b_j[j] = beta_t_1 * a_ij * dmvnorm(y_i[t+1, ], 
                                                       mean = m_list[[j]], 
                                                       sigma = cov_list[[j]])
                }
                
                beta_mat[t,l] = sum(b_j)
            }
        }
    }
    
    return(beta_mat)
}

log_likelihood_fnc <- function(par, par_index, y, id) {
    
    # Initial state probabilities
    init_logit = c(1, exp(pars[par_index$logit_pi][1]), 
                   exp(pars[par_index$logit_pi][2]))
    init = init_logit / sum(init_logit)
    
    m_1 = par[par_index$mu_1] # mean for state 1
    m_2 = par[par_index$mu_2] # mean for state 2
    m_3 = par[par_index$mu_3] # mean for state 3
    
    cov_1 = diag(exp(par[par_index$log_sig_1])) # covariance for state 1
    cov_2 = diag(exp(par[par_index$log_sig_2])) # covariance for state 2
    cov_3 = diag(exp(par[par_index$log_sig_3])) # covariance for state 3
    
    # Transition probability matrix
    q1   = exp(par[par_index$logit_t_p][1])
    q2   = exp(par[par_index$logit_t_p][2])
    q3   = exp(par[par_index$logit_t_p][3])
    q4   = exp(par[par_index$logit_t_p][4])
    q5   = exp(par[par_index$logit_t_p][5])
    q6   = exp(par[par_index$logit_t_p][6])
    
    Q = matrix(c(1,  q1, q2,
                 q3,   1, q4,
                 q5,  q6,  1), ncol=3, byrow=T)
    P = Q / rowSums(Q)
    
    # Parallelized computation of the log-likelihood
    log_total_val = foreach(i=unique(id), .combine='+') %dopar% {
        
        y_i = y[id == i, , drop = F]  # the matrix of responses for subject i
        
        d_1 = dmvnorm(y_i[1, ], mean = m_1, sigma = cov_1)
        d_2 = dmvnorm(y_i[1, ], mean = m_2, sigma = cov_2)
        d_3 = dmvnorm(y_i[1, ], mean = m_3, sigma = cov_3)
        
        f_i = init %*% diag(c(d_1, d_2, d_3))
        log_norm = 0
        val = 1
        
        for(t in 2:nrow(y_i)) {
            
            d_1 = dmvnorm(y_i[t, ], mean = m_1, sigma = cov_1)
            d_2 = dmvnorm(y_i[t, ], mean = m_2, sigma = cov_2)
            d_3 = dmvnorm(y_i[t, ], mean = m_3, sigma = cov_3)
            
            val = f_i %*% P %*% diag(c(d_1, d_2, d_3))
            
            norm_val = sqrt(sum(val^2))
            f_i = val / norm_val
            log_norm = log_norm + log(norm_val)
        }
        
        return(log(sum(f_i)) + log_norm)
    }
    
    return(log_total_val)
}
