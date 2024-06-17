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
            arma::vec norm_vec = dmvnorm(y_i.row(t), m_list(i), cov_list(i), false);
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

arma::mat forward_proc_it_c(arma::field<arma::vec> m_list, arma::field<arma::mat> cov_list,
                            arma::vec init, arma::mat P, arma::mat y_i) {
    
    arma::mat alpha_mat(y_i.n_rows, m_list.n_elem);
    
    for(int t = 0; t < y_i.n_rows; t++) {
        if(t == 0) {
            for(int l = 0; l < m_list.n_elem; l++) {
                arma::vec norm_vec = dmvnorm(y_i.row(t), m_list(l), cov_list(l), false);
                alpha_mat(t,l) = arma::as_scalar(init(l) * norm_vec);
            }
        } else {
            for(int l = 0; l < m_list.n_elem; l++) {
                arma::vec a_ji = P.col(l);
                arma::rowvec alpha_t_l = alpha_mat.row(t-1);
                
                double a_ji_alpha_prod = arma::as_scalar(alpha_t_l * a_ji);
                arma::vec norm_vec = dmvnorm(y_i.row(t), m_list(l), cov_list(l), false);
                
                alpha_mat(t,l) = a_ji_alpha_prod * arma::as_scalar(norm_vec);
            }
        }
    }
    
    return alpha_mat;
}

arma::mat backward_proc_it_c(arma::field<arma::vec> m_list, arma::field<arma::mat> cov_list,
                             arma::vec init, arma::mat P, arma::mat y_i) {
    
    arma::mat beta_mat(y_i.n_rows, m_list.n_elem);
    
    for(int t = y_i.n_rows - 1; t >= 0; t--) {
        if(t == y_i.n_rows - 1) {
            arma::rowvec beta_vec_1(beta_mat.n_cols, arma::fill::ones);
            beta_mat.row(t) = beta_vec_1;
        } else {
            for(int l = 0; l < m_list.n_elem; l++) {
                arma::vec b_j(m_list.n_elem, arma::fill::zeros);
                for(int j = 0; j < m_list.n_elem; j++) {
                    double beta_t_1 = beta_mat(t+1, j);
                    double a_ij = P(l,j);
                    
                    arma::vec norm_vec = dmvnorm(y_i.row(t+1), m_list(j), cov_list(j), false);
                    b_j(j) = beta_t_1 * a_ij * arma::as_scalar(norm_vec);
                }
                
                beta_mat(t,l) = arma::accu(b_j);
            }
        }
    }
    
    return beta_mat;
}

// [[Rcpp::export]]
arma::vec mu_s_update_c(const int s, arma::field<arma::mat> big_gamma, 
                      arma::mat y, arma::vec id) {
    
    arma::vec id_unique = arma::unique(id);
    arma::vec mu_hat_num(y.n_cols, arma::fill::zeros);
    double mu_hat_den = 0;

    for(int i = 0; i < id_unique.n_elem; i++) {
        arma::uvec sub_ind = arma::find(id == id_unique(i));
        arma::mat y_i = y.rows(sub_ind);
        arma::vec gamma_i_s = big_gamma(i).col(s-1);
        
        arma::mat diag_gamma = arma::diagmat(gamma_i_s); 
        
        arma::mat gamma_y = diag_gamma * y_i;
        
        arma::rowvec gamma_col_sum = arma::sum(gamma_y, 0);
        
        mu_hat_num += gamma_col_sum.t();
        mu_hat_den += arma::accu(gamma_i_s);
    }
    
    arma::vec mu_hat = mu_hat_num / mu_hat_den;
    
    return mu_hat;
}

// [[Rcpp::export]]
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
        arma::uvec sub_ind = arma::find(id == id_unique(i));
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

// [[Rcpp::export]]
arma::rowvec pi_s_update_c(arma::field<arma::mat> big_gamma, arma::vec id) {
    
    arma::rowvec pi_hat_num(3, arma::fill::zeros);
    double pi_hat_den = 0;
    
    arma::vec id_unique = arma::unique(id);
    
    for(int i = 0; i < id_unique.n_elem; i++) {
        arma::mat gamma_i = big_gamma(i);
        
        pi_hat_num = pi_hat_num + gamma_i.row(0);
        pi_hat_den += arma::accu(gamma_i.row(0));
    }
    
    arma::rowvec pi_hat = pi_hat_num / pi_hat_den;
    
    return pi_hat;
}

// [[Rcpp::export]]
arma::mat A_sm_update_c(arma::field<arma::field<arma::mat>> big_gamma,
                        arma::field<arma::field<arma::field<arma::field<arma::vec>>>> big_xi,
                        arma::field<arma::vec> id, int n_env) {
    
    arma::mat a_sm_num(3, 3, arma::fill::zeros);
    arma::mat a_sm_den(3, 3, arma::fill::zeros);
    
    arma::mat A_sm_hat(3, 3, arma::fill::zeros);
    arma::mat A_sm_hat2(3, 3, arma::fill::zeros);
    
    for(int ee = 0; ee < n_env; ee++) {
        arma::vec id_unique = arma::unique(id(ee));
        
        arma::mat a_sm_num2(3, 3, arma::fill::zeros);
        arma::mat a_sm_den2(3, 3, arma::fill::zeros);
        
        for(int i = 0; i < id_unique.n_elem; i++) {
            for(int s = 0; s < 3; s++) {
                for(int m = 0; m < 3; m++) {
                    arma::mat gamma_i = big_gamma(ee)(i);
                    arma::field<arma::field<arma::vec>> xi_i = big_xi(ee)(i);
                    
                    // subset the s column to be all BUT the last element
                    arma::vec gamma_i_col_s = gamma_i.col(s).subvec(0, gamma_i.n_rows - 2);
                    
                    a_sm_num(s,m) = a_sm_num(s,m) + arma::accu(xi_i(s)(m));
                    a_sm_den(s,m) = a_sm_den(s,m) + arma::accu(gamma_i_col_s);    
                    
                    a_sm_num2(s,m) = a_sm_num2(s,m) + arma::accu(xi_i(s)(m));
                    a_sm_den2(s,m) = a_sm_den2(s,m) + arma::accu(gamma_i_col_s);
                }
            }
        }
        
        A_sm_hat2 = A_sm_hat2 + (a_sm_num2 / a_sm_den2);
    }
    
    A_sm_hat = a_sm_num / a_sm_den;
    A_sm_hat2 = A_sm_hat2 / n_env;
    
    return A_sm_hat2;
    
}

// [[Rcpp::export]]
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
    
    arma::mat P = arma::reshape(prob_val, m_list.n_elem, m_list.n_elem);
    
    arma::vec init = par.elem(par_index(1) - 1);
    
    
    // After each calculation of omega, we can save these calculations for -----
    // updating the other parameters -------------------------------------------
    
    arma::field<arma::mat> big_gamma(id_unique.n_elem);
    arma::field<arma::field<arma::field<arma::vec>>> big_xi(id_unique.n_elem);
    
    for(int i = 0; i < id_unique.n_elem; i++) {
        arma::uvec sub_ind = arma::find(id == id_unique(i));
        arma::mat y_i = y.rows(sub_ind);
        int n_i = y_i.n_rows;
        
        // Calculate forward proc. & backward proc.
        arma::mat alpha_mat = forward_proc_it_c(m_list, cov_list, init, P, y_i);
        arma::mat beta_mat  = backward_proc_it_c(m_list, cov_list, init, P, y_i);
        
        // Calculate gamma for all time points
        arma::mat gamma_mat = gamma_calc_c(alpha_mat, beta_mat, m_list, cov_list, init, P, y_i);
        big_gamma(i) = gamma_mat;
        
        // pi calculation
        // Rcpp::Rcout << init.t() << std::endl;
        // Rcpp::Rcout << log(init.t()) << std::endl;
        pi_comp = pi_comp + arma::as_scalar(gamma_mat.row(0) * log(init));
        
        // transition prob calculation
        big_xi(i).set_size(m_list.n_elem);
        for(int l = 0; l < m_list.n_elem; l++) {
            big_xi(i)(l).set_size(m_list.n_elem);
            for(int j = 0; j < m_list.n_elem; j++) {
                arma::vec xi_time = xi_calc_c(l, j, alpha_mat, beta_mat, m_list, cov_list, init, P, y_i);
                big_xi(i)(l)(j) = xi_time;
                
                A_comp = A_comp + arma::accu(log(P(l,j)) * xi_time);
            }
        }
    
        // likelihood calculation
        for(int t = 0; t < n_i; t++) {
            for(int l = 0; l < m_list.n_elem; l++) {
                arma::vec norm_vec = dmvnorm(y_i.row(t), m_list(l), cov_list(l), true);
                
                like_comp += gamma_mat(t, l) * arma::as_scalar(norm_vec);
            }
        }
    }
    
    double Q_k = pi_comp + A_comp + like_comp;
    
    Rcpp::List omega_list = Rcpp::List::create(Q_k, big_gamma, big_xi);
    
    return omega_list;
}

// [[Rcpp::export]]
double log_likelihood_fnc_c(arma::vec par, arma::field<arma::uvec> par_index, 
                            arma::mat y, arma::vec id) {
    // par_index KEY: (0) t_p, (1) init, (2) mu_1, (3) mu_2, (4) mu_3, (5) Sig_1, 
    //                (6) Sig_2, (7) Sig_3
    
    // Initialize key pieces ---------------------------------------------------
    arma::vec init_t = par.elem(par_index(1) - 1);
    arma::rowvec init = init_t.t();
    
    arma::vec id_unique = arma::unique(id);
    
    arma::field<arma::vec> m_list(3);
    m_list(0) = par.elem(par_index(2) - 1);
    m_list(1) = par.elem(par_index(3) - 1);
    m_list(2) = par.elem(par_index(4) - 1);
    
    arma::field<arma::mat> cov_list(3);
    cov_list(0) = arma::reshape(par.elem(par_index(5) - 1), 5, 5);
    cov_list(1) = arma::reshape(par.elem(par_index(6) - 1), 5, 5);
    cov_list(2) = arma::reshape(par.elem(par_index(7) - 1), 5, 5);
    
    arma::mat P = arma::reshape(par.elem(par_index(0) - 1), 3, 3);
    //  ------------------------------------------------------------------------
    
    double log_total_val = 0;
    
    for(int i = 0; i < id_unique.n_elem; i++) {
        arma::uvec sub_ind = arma::find(id == id_unique(i));
        arma::mat y_i = y.rows(sub_ind);
        
        double d_1 = arma::as_scalar(dmvnorm(y_i.row(0), m_list(0), cov_list(0)));
        double d_2 = arma::as_scalar(dmvnorm(y_i.row(0), m_list(1), cov_list(1)));
        double d_3 = arma::as_scalar(dmvnorm(y_i.row(0), m_list(2), cov_list(2)));
        
        arma::vec d_diag = {d_1, d_2, d_3};
        arma::mat D = arma::diagmat(d_diag);
        
        arma::rowvec f_i = init * D;
        double log_norm = 0;
        arma::rowvec val(3);
        
        for(int t = 1; t < y_i.n_rows; t++) {
            d_1 = arma::as_scalar(dmvnorm(y_i.row(t), m_list(0), cov_list(0)));
            d_2 = arma::as_scalar(dmvnorm(y_i.row(t), m_list(1), cov_list(1)));
            d_3 = arma::as_scalar(dmvnorm(y_i.row(t), m_list(2), cov_list(2)));
            
            d_diag = {d_1, d_2, d_3};
            D = arma::diagmat(d_diag);
            
            val = f_i * P * D;
            double norm_val = arma::norm(val, 2);
            f_i = val / norm_val;
            log_norm = log_norm + log(norm_val);
        }
        
        log_total_val = log_total_val + log(arma::accu(f_i)) + log_norm;
    }

    return log_total_val;

}

// [[Rcpp::export]]
int test_fnc(arma::vec par, arma::field<arma::uvec> par_index, arma::mat y) {
    
    arma::mat M = {{1,2,3,4},
                    {5,6,7,8},
                    {9,10,11,12},
                    {13,14,15,16}};
    arma::vec M_s = M.col(1).subvec(0, M.n_rows - 2);
    // Rcpp::Rcout << M << std::endl;
    // Rcpp::Rcout << M_s << std::endl;
    
    arma::field<arma::mat> f(4);
    
    f(0) = arma::eye(3,3);
    f(1) = arma::eye(3,3);
    f(2) = arma::eye(3,3);
    f(3) = arma::eye(3,3);
    
    // Rcpp::Rcout << f.n_elem << std::endl;
    
    arma::mat l = {{1,2},
                   {3,4}};
    
    arma::mat k = {{5,6},
                   {7,8}};
    
    Rcpp::Rcout << l % k << std::endl;
    
    arma::mat alpha_beta_prod = l % k;
    arma::vec alpha_beta_sum  = arma::sum(alpha_beta_prod, 1);
    
    arma::vec alpha_l = l.col(0);
    arma::vec beta_l  = k.col(0);
    
    Rcpp::Rcout << alpha_l % beta_l << std::endl;
    Rcpp::Rcout << alpha_beta_sum << std::endl;
    Rcpp::Rcout << (alpha_l % beta_l) / alpha_beta_sum << std::endl;
    
    // arma::vec init = {1,2,3,4};
    // 
    // Rcpp::Rcout << arma::accu(log(M(0,1)) * init) << std::endl;
    // 
    // Rcpp::Rcout << arma::as_scalar(M.row(1) * log(init)) << std::endl;
    
    arma::field<arma::vec> m_list(3);
    m_list(0) = par.elem(par_index(2) - 1);
    m_list(1) = par.elem(par_index(3) - 1);
    m_list(2) = par.elem(par_index(4) - 1);
    
    arma::field<arma::mat> cov_list(3);
    cov_list(0) = arma::reshape(par.elem(par_index(5) - 1), 5, 5);
    cov_list(1) = arma::reshape(par.elem(par_index(6) - 1), 5, 5);
    cov_list(2) = arma::reshape(par.elem(par_index(7) - 1), 5, 5);
    
    // arma::vec d_1_vec = dmvnorm(y.row(0), m_list(0), cov_list(0));
    // double d_1 = arma::as_scalar(dmvnorm(y.row(0), m_list(0), cov_list(0)));
    // 
    // Rcpp::Rcout << d_1_vec << std::endl;
    // Rcpp::Rcout << d_1 << std::endl;
    // 
    // arma::rowvec t_1 = {3, 4};
    // Rcpp::Rcout << arma::norm(t_1, 2) << std::endl;
    // 
    // arma::vec init = par.elem(par_index(1) - 1); 
    // 
    // arma::mat g_temp = {{0.5,1,1.5},
    //                     {1,2,3}};
    // 
    // Rcpp::Rcout << g_temp.row(0) * log(init) << std::endl;
    // Rcpp::Rcout << arma::as_scalar(g_temp.row(0) * log(init)) << std::endl;
    
    return 0; 
}













