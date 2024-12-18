

#' infer_pseudotime
#'
#' A function to compute the mean of a vector
#' @param fit Expression matrix: cells * genes
#'
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }
#'
infer_pseudotime <- function(bt) {

  t <- posterior.mode(mcmc(rstan::extract(bt$fit, paste0("t"))[[paste0("t")]]))
  # t_samples <- rstan::extract(fit, "t")$t
  # posterior_mean <- mean(t_samples)  # 标量参数
  # posterior_mean
  bt$t <- t
  return(bt)
}


#' infer_mu0
#'
#' A function to compute the mean of a vector
#' @param bf Expression matrix: cells * genes
#'
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }
#'
infer_mu0 <- function(bt) {
  Y <- bt$y
  n_branch <- bt$n_branch
  para_mu0 <- c()
  for (i in 1:n_branch){
    p_mu0 <- posterior.mode(mcmc(rstan::extract(bt$fit)$all_mu0[,i,]))
    para_mu0 = rbind(para_mu0,p_mu0)
  }
  colnames(para_mu0) <- colnames(Y)
  rownames(para_mu0) <- paste0('branch_',1:n_branch)

  bt$para_mu0 <- para_mu0
  return(bt)
}

#' infer_k
#'
#' A function to compute the mean of a vector
#' @param bf Expression matrix: cells * genes
#'
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }
#'
infer_k <- function(bt) {
  Y <- bt$y
  n_branch <- bt$n_branch
  para_k <- c()
  for (i in 1:n_branch){
    p_k <- posterior.mode(mcmc(rstan::extract(bt$fit)$all_k[,i,]))
    para_k = rbind(para_k,p_k)
  }
  colnames(para_k) <- colnames(Y)
  rownames(para_k) <- paste0('branch_',1:n_branch)
  bt$para_k <- para_k
  return(bt)
}

#' infer_t0
#'
#' A function to compute the mean of a vector
#' @param bf Expression matrix: cells * genes
#'
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }
#'
infer_t0 <- function(bt) {
  Y <- bt$y
  n_branch <- bt$n_branch
  para_t0 <- c()
  for (i in 1:n_branch){
    p_t0 <- posterior.mode(mcmc(rstan::extract(bt$fit)$all_t0[,i,]))
    para_t0 = rbind(para_t0,p_t0)
  }
  colnames(para_t0) <- colnames(Y)
  rownames(para_t0) <- paste0('branch_',1:n_branch)

  bt$para_t0 <- para_t0
  return(bt)
}



#' branch_probability
#'
#' A function to compute the mean of a vector
#' @param bf Expression matrix: cells * genes
#'
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }
#'
#'
branch_probability <- function(bt){
  Y <- bt$y
  n_branch <- bt$n_branch
  pb <- txtProgressBar(min = 0, max = (nrow(Y)), style = 3)
  p_matrix <- matrix(nrow = 0,ncol = n_branch)
  for (i in 1:(nrow(Y))) {
    category_prob <- apply(matrix(rstan::extract(bt$fit)[['category_probs']][,i,]),2,mean)
    p_matrix <- rbind(p_matrix,category_prob)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  rownames(p_matrix) <- rownames(Y)
  colnames(p_matrix) <- paste0('p',1:n_branch)
  posterior_probs <- p_matrix

  bt$posterior_probs <- posterior_probs
  return(bt)
}



#' calculate_entropy
#'
#' A function to compute the mean of a vector
#' @param bf Expression matrix: cells * genes
#'
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }
#'
calculate_entropy <- function(bt){
  entropy <- apply(bt$posterior_probs,1,function(d){
    res <- 0
    for(i in 1:ncol(bt$posterior_probs))
    {
      d_num <- as.numeric(d[i])
      if(d_num!=0)
        res <- res + d_num*log(d_num,base = 2)
    }
    return (-res)
  } )
  bt$entropy <- entropy
  return(bt)
}


#' calculate_branches
#'
#' A function to compute the mean of a vector
#' @param bf Expression matrix: cells * genes
#' @param limit_prob 概率大于多少认为是这个分支的细胞
#'
#' @importFrom MCMCglmm posterior.mode
#' @importFrom rstan extract
#' @importFrom coda mcmc
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }
#'
calculate_branches <- function(bt,limit_prob = 0.9){
  branch_probs <- bt$posterior_probs
  n_branch <- bt$n_branch
  min_t0 <- min(bt$para_t0)

  all_branch_data <- list()
  branch_ind <- list()
  defined_ind <- c()
  for(i in 1:n_branch){
    branch_i <- which(branch_probs[,i] > limit_prob & bt$t > min_t0)
    branch_ind[[i]] <- branch_i
    defined_ind <- c(defined_ind,branch_i)
  }
  undefined_ind <- setdiff(1:bt$N,defined_ind)
  for (i in 1:n_branch) {
    comb <- list()
    comb$index <- c(undefined_ind,branch_ind[[i]])
    comb$t <- bt$t[comb$index]
    comb$expr <- bt$Y[comb$index,]
    all_branch_data[[i]] <- comb
  }
  bt$all_branch_data <- all_branch_data
  return(bt)
}

