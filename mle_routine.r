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

# -----------------------------------------------------------------------------
# THOUGHTS??
# I think we should be able to precompute alpha and beta for all time points, t,
# at the start of each EM step that way we aren't computing the recursive alg.
# for each successive time point 
# -----------------------------------------------------------------------------


forward_proc <- function(state_l, m_list, cov_list, init, P, y_i, t) {

    # m_list[[state_l]] = mean for state = l
    # cov_list[[state_l]] = covariance for state = l

    if(t == 1) {
        # Basecase
        alpha_l_1 = init[state_l] * dmvnorm(y_i[t, ],mean = m_list[[state_l]], 
                                                    sigma = cov_list[[state_l]])
        return(alpha_l_1)

    } else {
        # Recursive step


    }

}

backward_proc <- function(state_l, m_list, cov_list, init, P, y_i, t) {

    # Recursive procedure

    if(t == 1) {
        # Basecase

    } else {
        # Recursive step

    }


}