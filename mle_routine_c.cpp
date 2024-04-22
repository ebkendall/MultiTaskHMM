#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


arma::vec mu_s_update(const int s, arma::field<arma::mat> big_gamma, 
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

arma::mat Sigma_s_update(const int s, arma::field<arma::mat> big_gamma,
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

double pi_s_update(const int s, arma::field<arma::mat> big_gamma, arma::vec id) {
    
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

double A_sm_update(int ind_j, arma::field<arma::field<arma::mat>> big_gamma,
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
            
        }
    }
    
    return 0;
    
}

// A_sm_update <- function(ind_j, big_gamma, big_xi, id, n_env) {
//     
//     pos_matrix = matrix(c(0,1,2,3,0,4,5,6,0),nrow = 3, byrow = T)
//     ind_j_pos = which(pos_matrix == ind_j, arr.ind = T)
//     s = ind_j_pos[1,1]
//     m = ind_j_pos[1,2]
//     
//     a_sm_num <- 0
//     a_sm_den <- 0
//     
//     for(ee in 1:n_env) {
//         id_unique <- unique(id[[ee]])
//         for(i in 1:length(id_unique)) {
//             
//             gamma_i = big_gamma[[ee]][[i]]
//             xi_i = big_xi[[ee]][[i]]
//             
//             a_sm_num <- a_sm_num + sum(xi_i[[s]][[m]])
//             a_sm_den <- a_sm_den + sum(gamma_i[1:(nrow(gamma_i) - 1), s])
//             
//         }   
//     }
//     
//     A_sm_hat <- a_sm_num / a_sm_den
//     
//     return(A_sm_hat)
// }

















