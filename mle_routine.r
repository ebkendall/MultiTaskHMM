library(mvtnorm, quietly=T) 
library(foreach, quietly=T) 
library(doParallel, quietly=T)

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

baum_welch_one_environment <- function(par, par_index, y, id) {

    eps <- 1e-4

    # while (difference > eps) {
    # 
    # }
}


# -----------------------------------------------------------------------------
# THOUGHTS??
# I think we should be able to precompute alpha and beta for all time points, t,
# at the start of each EM step that way we aren't computing the recursive alg.
# for each successive time point 
# -----------------------------------------------------------------------------
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