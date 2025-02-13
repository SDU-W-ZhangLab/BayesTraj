# BayesTraj
BayesTraj is a probabilistic pseudotime framework. BayesTraj

- **A.** Inputs to BayesTraj

- **B.** Illustration of the Bayesian Model

- **C.** Trend of Marker Gene Changes During Differentiation
  
- **D.** Generation of Pseudo-time and Branch Probabilities

- **E.** Branch-specific gene detection

 ![fig1_画板 1](https://github.com/user-attachments/assets/48514e5e-3276-414d-9dd5-a176af41855c)




# Getting started

## Installation

```r
# install.packages("devtools")
devtools::install_github("SDU-W-ZhangLab/BayesTraj")
```

## Model fitting

To fit the pseudotimes, pass the appropriately prepared data (for example, **the normalized gene expression matrices `Y1` and `Y2`**) to the `BayesTraj` function:

```r
library(BayesTraj)
Y1 <- (branch_1_exp)
Y2 <- (branch_2_exp) # maker gene expression data 
bt <- BayesTraj(iter = 100, n_branch = 2,
                Y1 = Y1,
                Y2 = Y2,
                alpha = 0.1, eplison = 0.04,
                lower_k_s = 5, lower_k_t = 1, lower_t0 = 0.1, upper_b1_t0 = 0.8, lower_b2_t0 = 0.1,
                b1_mu0_s_mu = b1_mu0_s_mu, b1_mu0_s_sd = b1_mu0_s_sd,
                # b1_k_s_mu = b1_k_s_mu, b1_k_s_sd = b1_k_s_sd,
                b2_mu0_s_mu = b2_mu0_s_mu, b2_mu0_s_sd = b2_mu0_s_sd,
                b1_mu0_t_mu = b1_mu0_t_mu, b1_mu0_t_sd = b1_mu0_t_sd,
                b2_mu0_t_mu = b2_mu0_t_mu, b2_mu0_t_sd = b2_mu0_t_sd,
                lower_mu0_s = 0.2, upper_mu0_s = 0.5, lower_mu0_t = 0.05, upper_mu0_t = 0.4,
                chains = 1)
bt <- infer_pseudotime(bt)
bt <- infer_mu0(bt)
bt <- infer_k(bt)
bt <- infer_t0(bt)
bt <- branch_probability(bt)
bt <- calculate_branches(bt)
```

- **`infer_pseudotime`**  
  The `infer_pseudotime` function extracts the maximum-a-posteriori (MAP) estimates of the pseudotimes for each cell.

- **`infer_mu0`**  
  The `infer_mu0` function extracts the MAP estimates of the initial expression values (`mu_0`) across all genes.

- **`infer_k`**  
  The `infer_k` function extracts the MAP estimates of the slope parameters (`k`), indicating how gene expression changes over pseudotime.

- **`infer_t0`**  
  The `infer_t0` function extracts the MAP estimates of the breakpoint time (`t_0`), marking where expression patterns may shift.

- **`branch_probability`**  
  The `branch_probability` function calculates the posterior probability of branch membership for each cell.

- **`calculate_branches`**  
  The `calculate_branches` function applies these posterior estimates to finalize the branching structure and update the trajectory object accordingly.

For further usage options see the vignette. 
