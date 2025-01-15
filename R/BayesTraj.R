#' BayesTraj
#'
#' A function to compute the mean of a vector
#' @param Y1 description
#' @param Y2 description
#' @param Y3 description
#' @param Y4 description
#'
#' @param n_branch number of branches
#'
#' @param b1_mu0_s_mu description
#' @param b1_mu0_s_sd description
#' @param b2_mu0_s_mu description#'
#' @param b2_mu0_s_sd description
#' @param b3_mu0_s_mu description
#' @param b3_mu0_s_sd description
#' @param b4_mu0_s_mu description
#' @param b4_mu0_s_sd description
#' @param b1_mu0_t_mu description
#' @param b1_k_s_sd description
#' @param b1_k_t_mu description
#' @param b1_k_t_sd description
#' @param b2_k_s_mu description
#' @param b2_k_s_sd description
#' @param b2_k_t_mu description
#' @param b2_k_t_sd description
#' @param b3_k_s_mu description
#' @param b3_k_s_sd description
#' @param b3_k_t_mu description
#' @param b3_k_t_sd description
#' @param b4_k_s_mu description
#' @param b4_k_s_sd description
#' @param b4_k_t_mu description
#' @param b4_k_t_sd description
#' @param lower_k_s description
#' @param lower_k_t description
#' @param lower_t0 description
#' @param upper_t0 description
#' @param upper_b1_t0 description
#' @param lower_b2_t0 description
#' @param alpha description
#' @param eplison description
#'
#' @param lower_mu0_s description
#' @param lower_mu0_t description
#' @param upper_mu0_s description
#' @param upper_mu0_t description
#' @param cores escription
#'
#' @import rstan
#' @importFrom Rcpp loadModule
#' @importFrom stats coef lm prcomp
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }
#'

BayesTraj <- function(Y1=NULL,Y2=NULL,Y3=NULL,Y4=NULL,n_branch,
                                 b1_mu0_s_mu = NULL,b1_mu0_s_sd = NULL,b2_mu0_s_mu = NULL,b2_mu0_s_sd = NULL,
                                 b3_mu0_s_mu = NULL,b3_mu0_s_sd = NULL,b4_mu0_s_mu = NULL,b4_mu0_s_sd = NULL,
                                 b1_mu0_t_mu = NULL,b1_mu0_t_sd = NULL,b2_mu0_t_mu = NULL,b2_mu0_t_sd = NULL,
                                 b3_mu0_t_mu = NULL,b3_mu0_t_sd = NULL,b4_mu0_t_mu = NULL,b4_mu0_t_sd = NULL,
                                 b1_k_s_mu = NULL,b1_k_s_sd = NULL,b1_k_t_mu = NULL,b1_k_t_sd = NULL,
                                 b2_k_s_mu = NULL,b2_k_s_sd = NULL,b2_k_t_mu = NULL,b2_k_t_sd = NULL,
                                 b3_k_s_mu = NULL,b3_k_s_sd = NULL,b3_k_t_mu = NULL,b3_k_t_sd = NULL,
                                 b4_k_s_mu = NULL,b4_k_s_sd = NULL,b4_k_t_mu = NULL,b4_k_t_sd = NULL,
                                 lower_k_s = NULL,lower_k_t = NULL,lower_t0 = 0.2,upper_t0 = 0.8,
                                 upper_b1_t0 = 0.8,lower_b2_t0 = 0.2,
                                 alpha = 0.1,eplison = 0.04,
                                 lower_mu0_s = NULL, lower_mu0_t = NULL, upper_mu0_s = NULL,upper_mu0_t=NULL,
                                 cores = getOption("mc.cores",1L),
                                 ...) {
  library(rstan)

  Y <- cbind(Y1,Y2,Y3,Y4)
  Y <- Y[ , !duplicated(colnames(Y))]

  N <- nrow(Y1) # number of cells
  G <- ncol(Y)

  b1_G_s <- ncol(Y1)
  b1_G_t <- G - b1_G_s
  b2_G_s <- ifelse(is.null(Y2),0,ncol(Y2))
  b2_G_t <- G - b2_G_s
  b3_G_s <- ifelse(is.null(Y3),0,ncol(Y3))
  b3_G_t <- G - b3_G_s
  b4_G_s <- ifelse(is.null(Y4),0,ncol(Y4))
  b4_G_t <- G - b4_G_s

  if(is.null(b1_mu0_s_mu)) b1_mu0_s_mu <- as.array(rep(0.5, b1_G_s))
  if(is.null(b1_mu0_s_sd)) b1_mu0_s_sd <- as.array(rep(0.1, b1_G_s))
  if(is.null(b2_mu0_s_mu)) b2_mu0_s_mu <- as.array(rep(0.5, b2_G_s))
  if(is.null(b2_mu0_s_sd)) b2_mu0_s_sd <- as.array(rep(0.1, b2_G_s))
  if(is.null(b3_mu0_s_mu)) b3_mu0_s_mu <- as.array(rep(0.5, b3_G_s))
  if(is.null(b3_mu0_s_sd)) b3_mu0_s_sd <- as.array(rep(0.1, b3_G_s))
  if(is.null(b4_mu0_s_mu)) b4_mu0_s_mu <- as.array(rep(0.5, b4_G_s))
  if(is.null(b4_mu0_s_sd)) b4_mu0_s_sd <- as.array(rep(0.1, b4_G_s))

  if(is.null(b1_mu0_t_mu)) b1_mu0_t_mu <- as.array(rep(0.25, b1_G_t))
  if(is.null(b1_mu0_t_sd)) b1_mu0_t_sd <- as.array(rep(0.1, b1_G_t))
  if(is.null(b2_mu0_t_mu)) b2_mu0_t_mu <- as.array(rep(0.25, b2_G_t))
  if(is.null(b2_mu0_t_sd)) b2_mu0_t_sd <- as.array(rep(0.1, b2_G_t))
  if(is.null(b3_mu0_t_mu)) b3_mu0_t_mu <- as.array(rep(0.25, b3_G_t))
  if(is.null(b3_mu0_t_sd)) b3_mu0_t_sd <- as.array(rep(0.1, b3_G_t))
  if(is.null(b4_mu0_t_mu)) b4_mu0_t_mu <- as.array(rep(0.25, b4_G_t))
  if(is.null(b4_mu0_t_sd)) b4_mu0_t_sd <- as.array(rep(0.1, b4_G_t))

  if(is.null(b1_k_s_mu)) b1_k_s_mu <- as.array(rep(10, b1_G_s))
  if(is.null(b1_k_s_sd)) b1_k_s_sd <- as.array(rep(5, b1_G_s))
  if(is.null(b1_k_t_mu)) b1_k_t_mu <- as.array(rep(5, b1_G_t))
  if(is.null(b1_k_t_sd)) b1_k_t_sd <- as.array(rep(5, b1_G_t))

  if(is.null(b2_k_s_mu)) b2_k_s_mu <- as.array(rep(10, b2_G_s))
  if(is.null(b2_k_s_sd)) b2_k_s_sd <- as.array(rep(5, b2_G_s))
  if(is.null(b2_k_t_mu)) b2_k_t_mu <- as.array(rep(5, b2_G_t))
  if(is.null(b2_k_t_sd)) b2_k_t_sd <- as.array(rep(5, b2_G_t))

  if(is.null(b3_k_s_mu)) b3_k_s_mu <- as.array(rep(10, b3_G_s))
  if(is.null(b3_k_s_sd)) b3_k_s_sd <- as.array(rep(5, b3_G_s))
  if(is.null(b3_k_t_mu)) b3_k_t_mu <- as.array(rep(5, b3_G_t))
  if(is.null(b3_k_t_sd)) b3_k_t_sd <- as.array(rep(5, b3_G_t))

  if(is.null(b4_k_s_mu)) b4_k_s_mu <- as.array(rep(10, b4_G_s))
  if(is.null(b4_k_s_sd)) b4_k_s_sd <- as.array(rep(5, b4_G_s))
  if(is.null(b4_k_t_mu)) b4_k_t_mu <- as.array(rep(5, b4_G_t))
  if(is.null(b4_k_t_sd)) b4_k_t_sd <- as.array(rep(5, b4_G_t))

  if(is.null(lower_mu0_s)) lower_mu0_s <- 0.3
  if(is.null(upper_mu0_s)) upper_mu0_s <- 0.6

  if(is.null(lower_mu0_t)) lower_mu0_t <- 0.1
  if(is.null(upper_mu0_t)) upper_mu0_t <- 0.4

  if(is.null(lower_k_s)) lower_k_s <- 10
  if(is.null(lower_k_t)) lower_k_t <- 2

  lst1 <- rep(1,G)
  lst1[which(colnames(Y) %in% colnames(Y1))] <- 0
  lst2 <- rep(1,G)
  lst2[which(colnames(Y) %in% colnames(Y2))] <- 0
  lst3 <- rep(1,G)
  lst3[which(colnames(Y) %in% colnames(Y3))] <- 0
  lst4 <- rep(1,G)
  lst4[which(colnames(Y) %in% colnames(Y4))] <- 0

  data <- list(
    y = Y,
    lst1,lst2,lst3,lst4,
    b1_G_s = b1_G_s,b1_G_t = b1_G_t,b2_G_s=b2_G_s,b2_G_t=b2_G_t,b3_G_s=b3_G_s,b3_G_t=b3_G_t,b4_G_s=b4_G_s,b4_G_t=b4_G_t,
    lower_mu0_s = lower_mu0_s, lower_mu0_t = lower_mu0_t, upper_mu0_s = upper_mu0_s,upper_mu0_t=upper_mu0_t,
    lower_k_s = lower_k_s,lower_k_t = lower_k_t,lower_t0 = lower_t0,upper_t0 = upper_t0,
    upper_b1_t0 = upper_b1_t0,lower_b2_t0 = lower_b2_t0,
    G = G, N = N,n_branch = n_branch,
    alpha = alpha,eplison = eplison,
    b1_k_s_mu = as.array(b1_k_s_mu),b1_k_s_sd = as.array(b1_k_s_sd),b1_k_t_mu = as.array(b1_k_t_mu),b1_k_t_sd = as.array(b1_k_t_sd),
    b2_k_s_mu = as.array(b2_k_s_mu),b2_k_s_sd = as.array(b2_k_s_sd),b2_k_t_mu = as.array(b2_k_t_mu),b2_k_t_sd = as.array(b2_k_t_sd),
    b3_k_s_mu = as.array(b3_k_s_mu),b3_k_s_sd = as.array(b3_k_s_sd),b3_k_t_mu = as.array(b3_k_t_mu),b3_k_t_sd = as.array(b3_k_t_sd),
    b4_k_s_mu = as.array(b4_k_s_mu),b4_k_s_sd = as.array(b4_k_s_sd),b4_k_t_mu = as.array(b4_k_t_mu),b4_k_t_sd = as.array(b4_k_t_sd),
    b1_mu0_s_mu = as.array(b1_mu0_s_mu),b1_mu0_s_sd = as.array(b1_mu0_s_sd),b2_mu0_s_mu = as.array(b2_mu0_s_mu),b2_mu0_s_sd = as.array(b2_mu0_s_sd),
    b3_mu0_s_mu = as.array(b3_mu0_s_mu),b3_mu0_s_sd = as.array(b3_mu0_s_sd),b4_mu0_s_mu = as.array(b4_mu0_s_mu),b4_mu0_s_sd = as.array(b4_mu0_s_sd),
    b1_mu0_t_mu = as.array(b1_mu0_t_mu),b1_mu0_t_sd = as.array(b1_mu0_t_sd),b2_mu0_t_mu = as.array(b2_mu0_t_mu),b2_mu0_t_sd = as.array(b2_mu0_t_sd),
    b3_mu0_t_mu = as.array(b3_mu0_t_mu),b3_mu0_t_sd = as.array(b3_mu0_t_sd),b4_mu0_t_mu = as.array(b4_mu0_t_mu),b4_mu0_t_sd = as.array(b4_mu0_t_sd)
    )

  # stanfile <- "D:/study/2024/cellfate/stan/4_branch.stan"
  stanfile <- system.file("extdata", 'BayesTraj.stan', package = "BayesTraj")

  stanargs <- list(...)

  if (!("iter" %in% names(stanargs))) stanargs$iter <- 1000
  if (!("warmup" %in% names(stanargs))) stanargs$warmup <- stanargs$iter / 2
  if (!("chains" %in% names(stanargs))) stanargs$chains <- 1

  stanargs$data <- data
  stanargs$file <- stanfile
  stanargs$cores <- cores
  fit <- do.call(stan, stanargs)

  # 获取当前环境中的所有对象并打包为列表
  params <- mget(ls(environment()))
  BayesTraj <- structure(list(fit = fit,
                                    y = Y,
                                    params = params,
                                    lst1 = lst1,lst2=lst2,lst3=lst3,lst4=lst4,
                                    G = G, N = N,n_branch = n_branch,
                                    iter = stanargs$iter,
                                    chains = stanargs$chains,
                                    thin = stanargs$thin),
                               class = "BayesTraj")
  return(BayesTraj)
}

