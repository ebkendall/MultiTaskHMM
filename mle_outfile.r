# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)
library(latex2exp)

par_index = list(t_p = 1:9, init_pi = 10:12, 
                 mu_1 = 13:17, mu_2 = 18:22, mu_3 = 23:27, 
                 Sig_1 = 28:52, Sig_2 = 53:77, Sig_3 = 78:102)

shared_index = par_index$t_p

index_seeds = 1:100

load('Data/true_par.rda')

args = commandArgs(TRUE)
n_env = as.numeric(args[1])

labels <- c(TeX(r'($P(S1 \to S1)$)'), TeX(r'($P(S1 \to S2)$)'), TeX(r'($P(S1 \to S3)$)'), 
            TeX(r'($P(S2 \to S1)$)'), TeX(r'($P(S2 \to S2)$)'), TeX(r'($P(S2 \to S3)$)'),
            TeX(r'($P(S3 \to S1)$)'), TeX(r'($P(S3 \to S2)$)'), TeX(r'($P(S3 \to S3)$)'),
            TeX(r'(P(init S1))'), TeX(r'(P(init S2))'), TeX(r'(P(init S3))'),
            TeX(r'( (state 1) $\mu_1 [1]$)'), TeX(r'( (state 1) $\mu_1 [2]$)'), 
            TeX(r'( (state 1) $\mu_1 [3]$)'), TeX(r'( (state 1) $\mu_1 [4]$)'), 
            TeX(r'( (state 1) $\mu_1 [5]$)'),
            TeX(r'( (state 2) $\mu_2 [1]$)'), TeX(r'( (state 2) $\mu_2 [2]$)'), 
            TeX(r'( (state 2) $\mu_2 [3]$)'), TeX(r'( (state 2) $\mu_2 [4]$)'), 
            TeX(r'( (state 2) $\mu_2 [5]$)'),
            TeX(r'( (state 3) $\mu_3 [1]$)'), TeX(r'( (state 3) $\mu_3 [2]$)'), 
            TeX(r'( (state 3) $\mu_3 [3]$)'), TeX(r'( (state 3) $\mu_3 [4]$)'), 
            TeX(r'( (state 3) $\mu_3 [5]$)'),
            paste0("(state 1) Sigma[", 1:5, ", ", rep(1:5, each = 5), "]"),
            paste0("(state 2) Sigma[", 1:5, ", ", rep(1:5, each = 5), "]"),
            paste0("(state 3) Sigma[", 1:5, ", ", rep(1:5, each = 5), "]"))

par_est_mat = vector(mode = 'list', length = n_env)
for(e in 1:n_env) par_est_mat[[e]] = matrix(ncol = length(labels), nrow = length(index_seeds))


for(i in index_seeds) {
    load(paste0('Model_out/model_out_multi_e', n_env, '_ind', i, '.rda'))
    for(e in 1:n_env) {
        par_est_mat[[e]][i, ] = model_out[[1]][[e]]
    }   
}

# Plot and save the mcmc trace plots and histograms.
pdf(paste0('Plots/par_est_', n_env, '.pdf'))

for(e in 1:n_env) {
    VP <- vector(mode="list", length = length(labels) + 1)
    truth_par = true_par[[e]]
    for(r in 1:length(labels)) {
        # Adding the boxplots
        yVar = par_est_mat[[e]][,r]
        disc_type = rep(n_env, nrow(par_est_mat[[e]]))
        y_label = ""
        
        title_color = "black"
        if(r %in% shared_index) title_color = "red"
        
        plot_df = data.frame(yVar = yVar, disc_type = disc_type)
        VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar)) +
            geom_violin(trim=FALSE) +
            geom_boxplot(width=0.1) +
            ggtitle(labels[r]) +
            ylab(y_label) +
            xlab(paste0("Parameter Value: ", round(truth_par[r], 3))) +
            geom_hline(yintercept=truth_par[r], linetype="dashed", color = "red") +
            theme(text = element_text(size = 7),
                  plot.title = element_text(color=title_color, face="bold"))
    }
    
    for(j in 1:12) {
        if(j < 12) {
            if(j==1) {
                if(e == 1) {
                    # only print the shared indices once
                    grid.arrange(VP[[(j-1)*9 + 1]], VP[[(j-1)*9 + 2]], VP[[(j-1)*9 + 3]], 
                                 VP[[(j-1)*9 + 4]], VP[[(j-1)*9 + 5]], VP[[(j-1)*9 + 6]], 
                                 VP[[(j-1)*9 + 7]], VP[[(j-1)*9 + 8]], VP[[(j-1)*9 + 9]], ncol=3, nrow =3)       
                }
            } else {
                grid.arrange(VP[[(j-1)*9 + 1]], VP[[(j-1)*9 + 2]], VP[[(j-1)*9 + 3]], 
                             VP[[(j-1)*9 + 4]], VP[[(j-1)*9 + 5]], VP[[(j-1)*9 + 6]], 
                             VP[[(j-1)*9 + 7]], VP[[(j-1)*9 + 8]], VP[[(j-1)*9 + 9]], ncol=3, nrow =3)
            }
        } else {
            grid.arrange(VP[[100]], VP[[101]], VP[[102]], ncol=3, nrow =3)
        }
    }
}

dev.off()
