############################################################################################
# Package: JUMP
# Version: 1.0.1
# Data: 2023-05-22
# Authors: P. Lyu, Y. Li, X. Wen, and H. Cao.
############################################################################################
#' @title Replicability Analysis of High-Throughput Experiments
#'
#' @param pvals1 A numeric vector of p-values from study 1.
#' @param pvals2 A numeric vector of p-values from study 2.
#' @param alpha The FDR level to control, default is 0.05.
#' @param lambda The values of the tuning parameter to estimate pi_0. Must be in [0,1), default is seq(0.01, 0.8, 0.01).
#'
#' @return a list with the following elements:
#' \item{p.max}{The maximum of p-values across two studies.}
#' \item{jump.thr}{The estimated threshold of p.max to control FDR at level alpha.}
#'
#' @importFrom splines, Rcpp, RcppArmadillo, stats
#'
#' @export
#'
#' @examples
#' # Simulate p-values in two studies
#' m = 10000
#' h = sample(0:3, m, replace = TRUE, prob = c(0.9, 0.025, 0.025, 0.05))
#' states1 = rep(0, m); states2 = rep(0, m)
#' states1[which(h==2|h==3)] = 1; states2[which(h==1|h==3)] = 1
#' z1 = rnorm(m, states1*2, 1)
#' z2 = rnorm(m, states2*3, 1)
#' p1 = 1 - pnorm(z1); p2 = 1 - pnorm(z2)
#' # Run JUMP to identify replicable signals
#' res.jump = JUMP(p1, p2, alpha = 0.05)
#' sig.idx = which(res.jump$p.max <= res.jump$jump.thr)
#'
JUMP <- function(pvals1, pvals2, alpha = 0.05, lambda = seq(0.01, 0.8, 0.01)){
  m = length(pvals1)

  # Storey's method to estimate the proportions
  xi00.hat = c(); pi0.hat1 = c(); pi0.hat2 = c()
  for (i in 1:length(lambda)) {
    xi00.hat[i] <- sum(pvals1[1:m]>=lambda[i] & pvals2>=lambda[i]) / (m*(1-lambda[i])^2)
    pi0.hat1[i] <- sum(pvals1>=lambda[i]) / (m * (1-lambda[i]))
    pi0.hat2[i] <- sum(pvals2>=lambda[i]) / (m * (1-lambda[i]))
  }
  # fitting a cubic spline by Storey and Tibshirani (2003)
  fit1 <- lm(xi00.hat~., data.frame(cbind(xi00.hat, bs(lambda))))
  fit2 <- lm(pi0.hat1~., data.frame(cbind(pi0.hat1, bs(lambda))))
  fit3 <- lm(pi0.hat2~., data.frame(cbind(pi0.hat2, bs(lambda))))

  pred1 = predict(fit1)
  diff1 = abs(diff(pred1))

  pred2 = predict(fit2)
  diff2 = abs(diff(pred2))

  pred3 = predict(fit3)
  diff3 = abs(diff(pred3))

  xi00.hat = as.numeric(pred1[which.min(diff1)])
  xi01.hat = max(0, as.numeric(pred2[which.min(diff2)] - xi00.hat))
  xi10.hat = max(0, as.numeric(pred3[which.min(diff3)] - xi00.hat))

  xi.hat = c(xi00.hat, xi01.hat, xi10.hat)
  res <- jump_cutoff(pvals1, pvals2, xi.hat, alpha)

  return(list(p.max = res$p_max, jump.thr = res$thr_jump))
}


