library(mvtnorm, quietly=T) 
library(foreach, quietly=T) 
library(doParallel, quietly=T)

baum_welch_one_environment <- function(par, par_index, y, id) {

    eps <- 1e-4

    omega_k    <- omega_k_calc(par, par_index, y, id)
    omega_k_1 <- 0

    while (abs(omega_k - omega_k_1) > eps) {
        
        omega_k_1 = omega_k
        
        # Update pi
        for(s in 2:3) {
            par[par_index$init_pi][s-1] <- pi_s_update(s, par, par_index, y, id)
        }
        
        # Update A
        it_A <- 1
        for(s in 1:3) {
            for(m in 1:3) {
                if(s != m) {
                    par[par_index$t_p][it_A] <- A_sm_update(s, m, par, par_index, y, id)
                    it_A = it_A + 1
                }
            }
        }
        
        for(s in 1:3) {
            # Update mu
            par[par_index[[s + 2]]] = mu_s_update(s, par, par_index, y, id)
            
            # Update Sigma  
            par[par_index[[s + 5]]] = Sigma_s_update(s, par, par_index, y, id)
        }
        
        omega_k    <- omega_k_calc(par, par_index, y, id)
        
    }
}

mu_s_update <- function(s, par, par_index, y, id) {
    
    # Initialize key components -----------------------------------------------
    id_unique <- unique(id)
    
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
    # -------------------------------------------------------------------------
    
    mu_hat_num <- rep(0, ncol(y))
    mu_hat_den <- 0
    
    for(i in 1:length(id_unique)) {
        y_i = y[id == id_unique[i], ]
        n_i = nrow(y_i)
        
        for(t in 1:n_i) {
            gamma_s_i_t <- gamma_calc(t, s, m_list, cov_list, init, P, y_i)
            
            mu_hat_num <- mu_hat_num + gamma_s_i_t * y_i[t, ]
            mu_hat_den <- mu_hat_den + gamma_s_i_t
        }
    }
    
    mu_hat <- mu_hat_num / mu_hat_den
    
    return(mu_hat)
}

Sigma_s_update <- function(s, par, par_index, y, id) {
    
    # Initialize key components -----------------------------------------------
    id_unique <- unique(id)
    
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
    # -------------------------------------------------------------------------
    
    sigma_hat_num <- matrix(0, nrow = 5, ncol = 5)
    sigma_hat_den <- 0
    
    for(i in 1:length(id_unique)) {
        y_i = y[id == id_unique[i], ]
        n_i = nrow(y_i)
        
        for(t in 1:n_i) {
            gamma_s_i_t <- gamma_calc(t, s, m_list, cov_list, init, P, y_i)
            
            y_mu_diff <- matrix(c(y_i[t, ] - m_list[[s]]), ncol = 1)
            outer_prod <- y_mu_diff %*% t(y_mu_diff)
            
            sigma_hat_num <- sigma_hat_num + gamma_s_i_t * outer_prod
            sigma_hat_den <- sigma_hat_den + gamma_s_i_t
        }
    }
    
    Sigma_hat <- sigma_hat_num / sigma_hat_den
    
    return(Sigma_hat)
}

pi_s_update <- function(s, par, par_index, y, id) {
    
    # Initialize key components -----------------------------------------------
    id_unique <- unique(id)
    
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
    # -------------------------------------------------------------------------
    
    pi_hat_num <- 0
    pi_hat_den <- 0
    
    for(i in 1:length(id_unique)) {
        y_i = y[id == id_unique[i], ]
        
        gamma_i_sum <- rep(0, 3)
        for(l in 1:3) {
            gamma_i_sum[l] <- gamma_calc(1, l, m_list, cov_list, init, P, y_i)
        }
        
        pi_hat_num <- pi_hat_num + gamma_i_sum[s]
        pi_hat_den <- pi_hat_den + sum(gamma_i_sum)
    }
    
    pi_hat <- pi_hat_num / pi_hat_den
    
    return(pi_hat)
}

A_sm_update <- function(s, m, par, par_index, y, id) {
    
    # Initialize key components -----------------------------------------------
    id_unique <- unique(id)
    
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
    # -------------------------------------------------------------------------
    
    a_sm_num <- 0
    a_sm_den <- 0
    
    for(i in 1:length(id_unique)) {
        y_i = y[id == id_unique[i], ]
        n_i = nrow(y_i)
        
        for(t in 1:n_i) {
            if(t != 1) {
                a_sm_num <- a_sm_num + xi_calc(t, s, m, m_list, cov_list, init, P, y_i)
            }
            
            if(t != n_i) {
                a_sm_den <- a_sm_den + gamma_calc(t, s, m_list, cov_list, init, P, y_i)
            }
        }
    }
    
    A_sm_hat <- a_sm_num / a_sm_den
    
    return(A_sm_hat)
}

omega_k_calc <- function(par, par_index, y, id) {
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

    for(i in 1:length(id_unique)) {
        print(i)
        y_i = y[id == id_unique[i], ]
        n_i = nrow(y_i)

        # pi calculation
        print("pi calculator")
        for(l in 1:3) {
            pi_comp <- pi_comp + gamma_calc(1, l, m_list, cov_list, init, P, y_i) * log(init[l])
        }

        # transition prob calculation
        print("transition prob calculator")
        for(t in 2:n_i) {
            for(l in 1:3) {
                for(j in 1:3) {
                    # The diagonal components are functions of the others
                    if(j != l) {
                        print(paste0(l, ", ", j))
                        A_comp = A_comp + xi_calc(t, l, j, m_list, cov_list, init, P, y_i) * log(P[l, j])
                    }
                }
            }
        }

        # likelihood calculation
        print("likelihood component")
        for(t in 1:n_i) {
            for(l in 1:3) {
                print(paste0('t: ', t, ", l: ", l))
                like_comp = like_comp + gamma_calc(t, l, m_list, cov_list, init, P, y_i) * 
                                            dmvnorm(y_i[t, ], mean = m_list[[l]], sigma = cov_list[[l]], log = T)
            }
        }
    }

    Q_k = pi_comp + A_comp + like_comp

    return(Q_k)

}

gamma_calc <- function(t, l, m_list, cov_list, init, P, y_i) {
    alpha_t_vec = NULL
    beta_t_vec  = NULL

    for(i in 1:3) {
        alpha_t_vec[i] = forward_proc(t, i, m_list, cov_list, init, P, y_i)
        beta_t_vec[i]  = backward_proc(t, i, m_list, cov_list, init, P, y_i)
    }

    gamma_t_l = (alpha_t_vec[l] * beta_t_vec[l]) / sum(alpha_t_vec * beta_t_vec)
    
    return(gamma_t_l)
}

xi_calc <- function(t, l, j, m_list, cov_list, init, P, y_i) {
    alpha_t_vec = NULL
    beta_t1_vec = NULL
    b_t1_vec    = NULL

    for(i in 1:3) {
        alpha_t_vec[i] = forward_proc(t, i, m_list, cov_list, init, P, y_i)
        beta_t1_vec[i] = backward_proc(t+1, i, m_list, cov_list, init, P, y_i)
        b_t1_vec[i]    = dmvnorm(y_i[t+1, ], mean = m_list[[i]], sigma = cov_list[[i]])
    }

    xi_numerator   = alpha_t_vec[l] * P[l,j] * beta_t1_vec[j] * b_t1_vec[j]
    xi_denominator = 0
    for(k in 1:3) {
        for(w in 1:3) {
            xi_denominator = xi_denominator + alpha_t_vec[k] * P[k,w] * beta_t1_vec[w] * b_t1_vec[w]
        }
    }

    xi_t_l = xi_numerator / xi_denominator

    return(xi_t_l)
}

forward_proc <- function(t, l, m_list, cov_list, init, P, y_i) {

    # m_list[[l]] = mean for state = l
    # cov_list[[l]] = covariance for state = l

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
            alpha_sum = alpha_sum + P[l,j] * forward_proc(t-1, j, m_list, cov_list, init, P, y_i)
        }

        return(b_l_t * alpha_sum)
    }
}

backward_proc <- function(t, l, m_list, cov_list, init, P, y_i) {

    # m_list[[l]] = mean for state = l
    # cov_list[[l]] = covariance for state = l

    if(t == nrow(y_i)) {
        # Basecase
        return(1)

    } else {
        # Recursive step
        beta_sum = 0
        for(j in 1:3) {
            b_j_t = dmvnorm(y_i[t, ], mean = m_list[[j]], sigma = cov_list[[j]])

            beta_sum = beta_sum + P[l,j] * b_j_t * backward_proc(t+1, j, m_list, cov_list, init, P, y_i)
        }

        return(beta_sum)

    }
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