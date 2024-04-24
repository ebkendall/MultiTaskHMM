#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

arma::mat gamma_calc_c(arma::mat alpha_mat, arma::mat beta_mat, 
                       arma::field<arma::vec> m_list, arma::field<arma::mat> cov_list, 
                       arma::vec init, arma::mat P, arma::mat y_i) {
    
    arma::mat gamma_mat(y_i.n_rows, m_list.n_elem);
    
    arma::mat alpha_beta_prod = alpha_mat % beta_mat;
    arma::vec alpha_beta_sum  = arma::sum(alpha_beta_prod, 1);
    
    for(int l = 0; l < gamma_mat.n_cols; l++) {
        arma::vec alpha_l = alpha_mat.col(l);
        arma::vec beta_l  = beta_mat.col(l);
        
        gamma_mat.col(l) = (alpha_l % beta_l) / alpha_beta_sum;
    }
    
    return gamma_mat;
} 

arma::vec xi_calc_c(int l, int j, arma::mat alpha_mat, arma::mat beta_mat, 
                    arma::field<arma::vec> m_list, arma::field<arma::mat> cov_list, 
                    arma::vec init, arma::mat P, arma::mat y_i) {
    
    arma::vec xi_vec(y_i.n_rows - 1, arma::fill::zeros);
    
    for(int t = 1; t < y_i.n_rows; t++){
        
        arma::vec b_t(m_list.n_elem, arma::fill::zeros);
        for(int i = 0; i < m_list.n_elem; i++){
            arma::vec norm_vec = dmvnorm(y_i.row(t), m_list(i), cov_list(i));
            b_t(i) = arma::as_scalar(norm_vec);
        }
        
        double xi_numerator = alpha_mat(t-1,l)*P(l,j)*beta_mat(t, j)*b_t(j);
        double xi_denominator = 0;
        
        for(int k = 0; k < m_list.n_elem; k++) {
            for(int w = 0; w < m_list.n_elem; w++) {
                xi_denominator += alpha_mat(t-1,k)*P(k,w)*beta_mat(t,w)*b_t(w);
            }
        }
        
        xi_vec(t-1) = xi_numerator / xi_denominator;
    }
    
    return xi_vec;
}

// forward_proc_it <- function(m_list, cov_list, init, P, y_i) {
//     
//     # m_list[[l]] = mean for state l
//     # cov_list[[l]] = covariance for state l
//     
//     alpha_mat = matrix(nrow = nrow(y_i), ncol = length(m_list))
//     
//     for(t in 1:nrow(y_i)) {
//         if(t == 1) {
//             for(l in 1:length(m_list)) {
//                 alpha_mat[t,l] = init[l] * dmvnorm(y_i[t, ], mean = m_list[[l]],
//                                                    sigma = cov_list[[l]])
//             }
//         } else {
//             for(l in 1:length(m_list)) {
//                 a_ji = c(P[,l])
//                 alpha_t_1 = alpha_mat[t-1, ]
//                 
//                 alpha_mat[t,l] = dmvnorm(y_i[t, ], mean = m_list[[l]], 
//                                          sigma = cov_list[[l]]) * sum(a_ji * alpha_t_1)
//             }
//         }
//     }
//     return(alpha_mat)
// }




arma::vec mu_s_update_c(const int s, arma::field<arma::mat> big_gamma, 
                      arma::mat y, arma::vec id) {
    
    arma::vec id_unique = arma::unique(id);
    arma::vec mu_hat_num(y.n_cols, arma::fill::zeros);
    double mu_hat_den = 0;

    for(int i = 0; i < id_unique.n_elem; i++) {
        arma::uvec sub_ind = arma::find(id == i);
        arma::mat y_i = y.rows(sub_ind);
        arma::vec gamma_i_s = big_gamma(i).col(s-1);
        
        arma::mat diag_gamma = arma::diagmat(gamma_i_s); 
        
        arma::mat gamma_y = diag_gamma * y_i;
        
        mu_hat_num += arma::sum(gamma_y, 0);
        mu_hat_den += arma::accu(gamma_i_s);
    }
    
    arma::vec mu_hat = mu_hat_num / mu_hat_den;
    
    return mu_hat;
}

arma::mat Sigma_s_update_c(const int s, arma::field<arma::mat> big_gamma,
                         arma::vec par, arma::field<arma::uvec> par_index,
                         arma::mat y, arma::vec id) {
    // par_index KEY: (0) t_p, (1) init, (2) mu_1, (3) mu_2, (4) mu_3, (5) Sig_1, 
    //                (6) Sig_2, (7) Sig_3
    
    // Initialize key pieces ---------------------------------------------------
    arma::vec id_unique = arma::unique(id);
    arma::field<arma::vec> m_list(3);
    m_list(0) = par.elem(par_index(2) - 1);
    m_list(1) = par.elem(par_index(3) - 1);
    m_list(2) = par.elem(par_index(4) - 1);
    //  ------------------------------------------------------------------------
    
    arma::mat sigma_hat_num(5, 5, arma::fill::zeros);
    double sigma_hat_den = 0;
    
    for(int i = 0; i < id_unique.n_elem; i++) {
        arma::uvec sub_ind = arma::find(id == i);
        arma::mat y_i = y.rows(sub_ind);
        int n_i = y_i.n_rows;
        
        for(int t = 0; t < n_i; t++) {
            double gamma_s_i_t = big_gamma(i)(t, s-1);
            
            arma::vec y_mu_diff = y_i.row(t).t() - m_list(s-1);
            arma::mat outer_prod = y_mu_diff * y_mu_diff.t();
            
            sigma_hat_num += gamma_s_i_t * outer_prod;
            sigma_hat_den += gamma_s_i_t;
        }
    }
    
    arma::mat Sigma_hat = sigma_hat_num / sigma_hat_den;
    
    return Sigma_hat;
}

double pi_s_update_c(const int s, arma::field<arma::mat> big_gamma, arma::vec id) {
    
    double pi_hat_num = 0;
    double pi_hat_den = 0;
    
    arma::vec id_unique = arma::unique(id);
    
    for(int i = 0; i < id_unique.n_elem; i++) {
        arma::mat gamma_i = big_gamma(i);
        
        pi_hat_num += gamma_i(0, s-1);
        pi_hat_den += arma::accu(gamma_i.row(0));
    }
    
    double pi_hat = pi_hat_num / pi_hat_den;
    
    return pi_hat;
}

double A_sm_update_c(int ind_j, arma::field<arma::field<arma::mat>> big_gamma,
                   arma::field<arma::field<arma::field<arma::field<arma::vec>>>> big_xi,
                   arma::field<arma::vec> id, int n_env) {
    // arma::mat pos_matrix = {{0, 1, 2},
    //                         {3, 0, 4}, 
    //                         {5, 6, 0}};
    arma::mat pos_mat = {{0, 0, 1, 1, 2, 2},
                         {1, 2, 0, 2, 0, 1}}; 
    arma::vec ind_j_pos = pos_mat.col(ind_j-1);
    int s = ind_j_pos(0);
    int m = ind_j_pos(1);
    
    double a_sm_num = 0;
    double a_sm_den = 0;
    
    for(int ee = 0; ee < n_env; ee++) {
        arma::vec id_unique = arma::unique(id(ee));
        for(int i = 0; i < id_unique.n_elem; i++) {
            arma::mat gamma_i = big_gamma(ee)(i);
            arma::field<arma::field<arma::vec>> xi_i = big_xi(ee)(i);
            
            // subset the s column to be all BUT the last element
            arma::vec gamma_i_col_s = gamma_i.col(s).subvec(0, gamma_i.n_rows - 2);
            
            a_sm_num += arma::accu(xi_i(s)(m));
            a_sm_den += arma::accu(gamma_i_col_s);
        }
    }
    
    double A_sm_hat = a_sm_num / a_sm_den;
    
    return A_sm_hat;
    
}

Rcpp::List omega_k_calc_c(arma::vec par, arma::field<arma::uvec> par_index, 
                      arma::mat y, arma::vec id) {
    // par_index KEY: (0) t_p, (1) init, (2) mu_1, (3) mu_2, (4) mu_3, (5) Sig_1, 
    //                (6) Sig_2, (7) Sig_3
    
    // Initialize key pieces ---------------------------------------------------
    arma::vec id_unique = arma::unique(id);
    arma::field<arma::vec> m_list(3);
    m_list(0) = par.elem(par_index(2) - 1);
    m_list(1) = par.elem(par_index(3) - 1);
    m_list(2) = par.elem(par_index(4) - 1);
    //  ------------------------------------------------------------------------
    
    double pi_comp = 0;
    double A_comp = 0;
    double like_comp = 0;
    
    arma::field<arma::mat> cov_list(3);
    cov_list(0) = arma::reshape(par.elem(par_index(5) - 1), 5, 5);
    cov_list(1) = arma::reshape(par.elem(par_index(6) - 1), 5, 5);
    cov_list(2) = arma::reshape(par.elem(par_index(7) - 1), 5, 5);
    
    arma::vec prob_val = par.elem(par_index(0) - 1);
    
    arma::mat P = {{1 - prob_val(0) - prob_val(1), prob_val(0), prob_val(1)},
                   {prob_val(2), 1 - prob_val(2) - prob_val(3), prob_val(3)},
                   {prob_val(4), prob_val(5), 1 - prob_val(4) - prob_val(5)}};
    
    arma::vec init_val = par.elem(par_index(1) - 1);
    arma::vec init = {1-arma::accu(init_val), init_val(0), init_val(1)};
    
    // After each calculation of omega, we can save these calculations for -----
    // updating the other parameters -------------------------------------------
    
    arma::field<arma::vec> big_gamma(id_unique.n_elem);
    arma::field<arma::field<arma::vec>> big_xi(id_unique.n_elem);
    
    for(int i = 0; i < id_unique.n_elem; i++) {
        arma::uvec sub_ind = arma::find(id == i);
        arma::mat y_i = y.rows(sub_ind);
        int n_i = y_i.n_rows;
        
        // Calculate forward proc. & backward proc.
        
        // Calculate gamma for all time points
        
        // pi calculation
        
        // transition prob calculation
        
        // likelihood calculation
    }
    
    double Q_k = pi_comp + A_comp + like_comp;
    
    Rcpp::List omega_list = Rcpp::List::create(Q_k, big_gamma, big_xi);
    
    return omega_list;
}



// [[Rcpp::export]]
int test_fnc() {
    
    arma::mat M = {{1,2,3,4},
                    {5,6,7,8},
                    {9,10,11,12},
                    {13,14,15,16}};
    arma::vec M_s = M.col(1).subvec(0, M.n_rows - 2);
    Rcpp::Rcout << M << std::endl;
    Rcpp::Rcout << M_s << std::endl;
    
    arma::field<arma::mat> f(4);
    
    f(0) = arma::eye(3,3);
    f(1) = arma::eye(3,3);
    f(2) = arma::eye(3,3);
    f(3) = arma::eye(3,3);
    
    Rcpp::Rcout << f.n_elem << std::endl;
    
    arma::mat l = {{1,2},
                   {3,4}};
    
    arma::mat k = {{5,6},
                   {7,8}};
    
    Rcpp::Rcout << l % k << std::endl;
    
    return 0; 
}













