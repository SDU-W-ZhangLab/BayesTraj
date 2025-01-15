#' Mean
#'
#' A function to compute the mean of a vector
#' @param x A numeric vector
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }
#'

BayesFates1 <- function(Y1=NULL,Y2=NULL,Y3=NULL,Y4=NULL,
                       b1_mu0_s_mu = NULL,b1_mu0_s_sd = NULL,b2_mu0_s_mu = NULL,b2_mu0_s_sd = NULL,
                       b3_mu0_s_mu = NULL,b3_mu0_s_sd = NULL,b4_mu0_s_mu = NULL,b4_mu0_s_sd = NULL,
                       b1_mu0_t_mu = NULL,b1_mu0_t_sd = NULL,b2_mu0_t_mu = NULL,b2_mu0_t_sd = NULL,
                       b3_mu0_t_mu = NULL,b3_mu0_t_sd = NULL,b4_mu0_t_mu = NULL,b4_mu0_t_sd = NULL,
                       b1_k_s_mu = NULL,b1_k_s_sd = NULL,b1_k_t_mu = NULL,b1_k_t_sd = NULL,
                       b2_k_s_mu = NULL,b2_k_s_sd = NULL,b2_k_t_mu = NULL,b2_k_t_sd = NULL,
                       b3_k_s_mu = NULL,b3_k_s_sd = NULL,b3_k_t_mu = NULL,b3_k_t_sd = NULL,
                       b4_k_s_mu = NULL,b4_k_s_sd = NULL,b4_k_t_mu = NULL,b4_k_t_sd = NULL,
                       lower_k_s = NULL,lower_k_t = NULL,
                       lower_t0 = 0.2,upper_t0 = 0.8,
                       upper_b1_t0 = 0.8,lower_b2_t0 = 0.2,
                       alpha = 0.1,eplison = 0.04,
                       lower_mu0_s = NULL, lower_mu0_t = NULL, upper_mu0_s = NULL,upper_mu0_t=NULL,
                       student_df = 10,n_branch,gene_pairs = list(),cores = getOption("mc.cores",1L),
                       ...) {
  library(rstan)

  # 将所有基因合并
  Y <- cbind(Y1,Y2,Y3,Y4,Y5,Y6)
  Y <- Y[ , !duplicated(colnames(Y))]
  ## Now sanitize the input
  N <- nrow(Y1) # number of cells 细胞数量
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
  if(is.null(b1_mu0_s_sd)) b1_mu0_s_sd <- as.array(rep(0.01, b1_G_s))
  if(is.null(b2_mu0_s_mu)) b2_mu0_s_mu <- as.array(rep(0.5, b2_G_s))
  if(is.null(b2_mu0_s_sd)) b2_mu0_s_sd <- as.array(rep(0.01, b2_G_s))
  if(is.null(b3_mu0_s_mu)) b3_mu0_s_mu <- as.array(rep(0.5, b3_G_s))
  if(is.null(b3_mu0_s_sd)) b3_mu0_s_sd <- as.array(rep(0.01, b3_G_s))
  if(is.null(b4_mu0_s_mu)) b4_mu0_s_mu <- as.array(rep(0.5, b4_G_s))
  if(is.null(b4_mu0_s_sd)) b4_mu0_s_sd <- as.array(rep(0.01, b4_G_s))

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

  if(is.null(lower_mu0_s)) lower_mu0_s <- 0.4
  if(is.null(upper_mu0_s)) upper_mu0_s <- 0.6

  if(is.null(lower_mu0_t)) lower_mu0_t <- 0.1
  if(is.null(upper_mu0_t)) upper_mu0_t <- 0.4

  if(is.null(lower_k_s)) lower_k_s <- 10
  if(is.null(lower_k_t)) lower_k_t <- 2

  G_d <- ifelse(is.null(Yd),0,ncol(Yd))
  if(is.null(Yd)){
    Yd <- matrix(NA,nrow = N,ncol = 0)
  }

  # if(is.null(b1_mu0_d_mu)) b1_mu0_d_mu <- as.array(rep(0.4, G_d))
  # if(is.null(b1_mu0_d_sd)) b1_mu0_d_sd <- as.array(rep(0.05, G_d))
  # if(is.null(b2_mu0_d_mu)) b2_mu0_d_mu <- as.array(rep(0.4, G_d))
  # if(is.null(b2_mu0_d_sd)) b2_mu0_d_sd <- as.array(rep(0.05, G_d))
  # if(is.null(b3_mu0_d_mu)) b3_mu0_d_mu <- as.array(rep(0.4, G_d))
  # if(is.null(b3_mu0_d_sd)) b3_mu0_d_sd <- as.array(rep(0.05, G_d))
  # if(is.null(b4_mu0_d_mu)) b4_mu0_d_mu <- as.array(rep(0.4, G_d))
  # if(is.null(b4_mu0_d_sd)) b4_mu0_d_sd <- as.array(rep(0.05, G_d))
  #
  # if(is.null(b1_k_d_mu)) b1_k_d_mu <- as.array(rep(10, G_d))
  # if(is.null(b1_k_d_sd)) b1_k_d_sd <- as.array(rep(10, G_d))
  # if(is.null(b2_k_d_mu)) b2_k_d_mu <- as.array(rep(10, G_d))
  # if(is.null(b2_k_d_sd)) b2_k_d_sd <- as.array(rep(10, G_d))
  # if(is.null(b3_k_d_mu)) b3_k_d_mu <- as.array(rep(10, G_d))
  # if(is.null(b3_k_d_sd)) b3_k_d_sd <- as.array(rep(10, G_d))
  # if(is.null(b4_k_d_mu)) b4_k_d_mu <- as.array(rep(10, G_d))
  # if(is.null(b4_k_d_sd)) b4_k_d_sd <- as.array(rep(10, G_d))

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
    lower_k_s = lower_k_s,lower_k_t = lower_k_t,
    lower_t0 = lower_t0,upper_t0 = upper_t0,
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
    b3_mu0_t_mu = as.array(b3_mu0_t_mu),b3_mu0_t_sd = as.array(b3_mu0_t_sd),b4_mu0_t_mu = as.array(b4_mu0_t_mu),b4_mu0_t_sd = as.array(b4_mu0_t_sd),
    student_df = student_df)

  # stanfile <- "D:/study/2024/cellfate/stan/4_branch.stan"
  stanfile <- system.file("extdata", 'BayesFates.stan', package = "mypackage")

  stanargs <- list(...)

  stanargs$data <- data
  stanargs$file <- stanfile
  stanargs$cores <- cores
  fit <- do.call(stan, stanargs)

  # 获取当前环境中的所有对象并打包为列表
  params <- mget(ls(environment()))
  multiple_branch <- structure(list(fit = fit,
                                    y = Y,
                                    params = params,
                                    G = G, N = N,n_branch = n_branch,
                                    iter = stanargs$iter,
                                    chains = stanargs$chains,
                                    thin = stanargs$thin),
                               class = "multiple_branch")
  return(multiple_branch)
}
