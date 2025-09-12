
#' Get the covariance estimates for the global factors
#'
#' @description
#' This function generates the covariance estimates for the global factors
#' at time \eqn{t}.
#'
#' @param object An S3 object of class 'multi_result' created by [multilevel()].
#' @param t An integer specifying the time
#'
#' @return An \eqn{r_{0} \times r_{0}} covariance matrix.
#' @import Matrix
#' @export
#'
#' @examples
#'
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' vcov <- vcov_global_factor(est_multi, t = est_multi$T / 2)
vcov_global_factor <- function(object, t) {

  G <- object$G
  F <- object$F
  Gamma <- object$Gamma
  Lambda <- object$Lambda
  R <- object$R
  N <- object$N
  Ni <- object$Ni
  T <- object$T
  r0 <- object$r0
  ri <- object$ri
  Resid <- object$Resid

  Var_G <- matrix(0, r0, r0)
  mI <- replicate(n = R, expr = diag(1, r0, r0), simplify = FALSE)
  mI <- do.call(rbind, mI)

  Ci <- list()
  Theta_list <- list()

  for (i in 1:R) {
    Ii <- rbind(diag(1, r0, r0), matrix(0, ri[i], r0))
    if (ri[i] > 0) {
      Thetai <- cbind(Gamma[[i]], Lambda[[i]])
    } else {
      Thetai <- Gamma[[i]]
    }
    Theta_list[[i]] <- Thetai / sqrt(Ni[i])
    Ci[[i]] <- sqrt(N / Ni[i]) * t(Ii) %*% solve(t(Thetai) %*% Thetai / Ni[i])
  }
  Theta_bd <- Matrix::bdiag(Theta_list)
  C_bd <- Matrix::bdiag(Ci)
  Sigma_e_thres <- vcov_POET(object)
  D1 <- t(Theta_bd) %*% Sigma_e_thres %*% Theta_bd
  Var_G <- (1 / R^2) * t(mI) %*% C_bd %*% D1 %*% t(C_bd) %*% mI
  return(as.matrix(Var_G))
}

#' Get the covariance estimates for the global factor loadings
#'
#' @description This function generates the covariance estimates
#' for the global factor loadings for the \eqn{j}-th individual in block \eqn{i}.
#'
#' @param object An S3 object of class 'multi_result' created by [multilevel()].
#' @param i An integer indicating the \eqn{i}-th block.
#' @param j An integer indicating the \eqn{j}-th individual in the \eqn{i}-th block.
#'
#' @return An \eqn{r_{0} \times r_{0}} covariance matrix.
#'
#' @export
#'
#' @examples
#'
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' vcov_gamma_11 <- vcov_global_loading(est_multi, i = 1, j = 1)
vcov_global_loading <- function(object, i, j) {
  G <- object$G
  F <- object$F
  Gamma <- object$Gamma
  Lambda <- object$Lambda
  R <- object$R
  N <- object$N
  Ni <- object$Ni
  T <- object$T
  r0 <- object$r0
  ri <- object$ri
  Resid <- object$Resid
  block <- i
  ind <- j

  Xij <- matrix(0, T, r0)
  if (ri[block] > 0) {
    lambda_ij <- matrix(Lambda[[block]][ind, ], ncol = 1)
  }
  gamma_ij <- matrix(Gamma[[block]][ind, ], ncol = 1)
  Ei <- Resid[[block]]
  for (s in 1:T) {
    if (ri[i] > 0) {
      Xij[s, ] <- G[s, ] * c(t(lambda_ij) %*% matrix(F[[block]][s, ], ncol = 1) + Ei[s, ind])
    } else {
      Xij[s, ] <- G[s, ] * c(Ei[s, ind])
    }
  }
  fit <- stats::lm(Xij ~ 1)
  Var_gammaij <- sandwich::NeweyWest(fit, sandwich = FALSE)

  return(Var_gammaij)
}

#' Get the covariance estimates for the local factors
#'
#' @description This function generates the covariance estimates
#' for the local factors in block \eqn{i} at time \eqn{t}.
#'
#' @param object An S3 object of class 'multi_result' created by multilevel().
#' @param i An integer indicating the \eqn{i}-th block.
#' @param t An integer specifying the time point.
#'
#' @return An \eqn{r_{i} \times r_{i}} covariance matrix.
#' @import Matrix
#' @export
#'
#' @examples
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' vcov_local_factor_11 <- vcov_local_factor(est_multi, i = 1, t = 1)
vcov_local_factor <- function(object, i, t) {
  G <- object$G
  F <- object$F
  Gamma <- object$Gamma
  Lambda <- object$Lambda
  R <- object$R
  N <- object$N
  Ni <- object$Ni
  T <- object$T
  r0 <- object$r0
  ri <- object$ri
  Resid <- object$Resid
  Y_block <- object$Y_list
  block <- i

  IIi <- rbind(matrix(0, r0, ri[block]), diag(1, ri[block], ri[block]))
  Yi_hat <- Y_block[[block]] - G %*% t(Gamma[[block]])
  eig <- eigen(Yi_hat %*% t(Yi_hat) / (Ni[block] * T))
  Upsi_hat <- matrix(0, ri[block], ri[block])
  diag(Upsi_hat) <- eig$values[1:ri[block]]

  Ci <- list()
  Theta_list <- list()
  Ci_dot <- list()
  for (i in 1:R) {
    Ii <- rbind(diag(1, r0, r0), matrix(0, ri[i], r0))
    if (ri[i] > 0) {
      Thetai <- cbind(Gamma[[i]], Lambda[[i]])
    } else {
      Thetai <- Gamma[[i]]
    }
    Theta_list[[i]] <- Thetai / sqrt(Ni[i])
    Ci[[i]] <- sqrt(N / Ni[i]) * t(Ii) %*% solve(t(Thetai) %*% Thetai / Ni[i])
    if (i != block) {
      Ci_dot[[i]] <- matrix(0, ri[block], r0 + ri[i])
    } else {
      Ci_dot[[i]] <- t(IIi)
    }
  }
  Theta_bd <- Matrix::bdiag(Theta_list)
  C_bd <- Matrix::bdiag(Ci)
  Sigma_e_thres <- vcov_POET(object)
  D1 <- t(Theta_bd) %*% Sigma_e_thres %*% Theta_bd
  mI <- replicate(n = R, expr = diag(1, r0, r0), simplify = FALSE)
  mI <- do.call(rbind, mI)
  Ci_dot <- do.call(cbind, Ci_dot) - sqrt(Ni[block] / N) * (1 / R) *
    (t(Lambda[[block]]) %*% Gamma[[block]] / Ni[block]) %*% t(mI) %*% C_bd
  Var_Fi <- solve(Upsi_hat) %*% Ci_dot %*% D1 %*% t(Ci_dot) %*% solve(Upsi_hat)

  return(as.matrix(Var_Fi))
}

#' Get the covariance estimates for the local factor loadings
#'
#' @description This function generates the covariance estimates
#' for the local loadings for the \eqn{j}-th individual in block \eqn{i}.
#'
#' @param object An S3 object of class 'multi_result' created by multilevel().
#' @param i An integer indicating the \eqn{i}-th block.
#' @param j An integer indicating the \eqn{j}-th individual in the \eqn{i}-th block.
#'
#' @return An \eqn{r_{i} \times r_{i}} covariance matrix.
#'
#' @export
#'
#' @examples
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' vcov_local_loading_11 <- vcov_local_loading(est_multi, i = 1, j = 1)
vcov_local_loading <- function(object, i, j) {

  F <- object$F
  Lambda <- object$Lambda
  T <- object$T
  Resid <- object$Resid
  Lambda_ij_hat <- Lambda[[i]][j, ]
  FEij <- F[[i]] * c(Resid[[i]][, j])
  fit <- stats::lm(FEij ~ 1)

  Var_FEij <- sandwich::NeweyWest(fit, sandwich = FALSE)

  return(Var_FEij)
}

#' Get the variance estimates of the global component
#'
#' @description This function generates the variance estimates of the
#' global component for the \eqn{j}-th individual in block \eqn{i} at time \eqn{t}.
#'
#' @param object An S3 object of class 'multi_result' created by multilevel().
#' @param i An integer indicating the \eqn{i}-th block.
#' @param j An integer indicating the \eqn{j}-th individual in the \eqn{i}-th block.
#' @param t An integer indicating the time.
#'
#' @return The variance of the global component.
#' @import Matrix
#' @export
#'
#' @examples
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' vcov_global_comp_ijt <- vcov_global_comp(est_multi, i = 1, j = 1, t = 1)
vcov_global_comp <- function(object, i, j, t) {
  G <- object$G
  F <- object$F
  Gamma <- object$Gamma
  Lambda <- object$Lambda
  R <- object$R
  N <- object$N
  Ni <- object$Ni
  T <- object$T
  r0 <- object$r0
  ri <- object$ri
  Resid <- object$Resid
  block <- i
  ind <- j

  D1t <- vcov_global_factor(object, t)
  D2ij <- vcov_global_loading(object, i = block, j = ind)

  gamma_ij <- matrix(Gamma[[block]][ind, ], ncol = 1)
  Gt <- matrix(G[t, ], ncol = 1)

  Var_gc_ijt <- (t(gamma_ij) %*% D1t %*% gamma_ij / N + t(Gt) %*% D2ij %*% Gt / T) * min(min(Ni),T)

  return(c(Var_gc_ijt))
}

#' Get the variance estimates of the local component
#'
#' @description This function generates the variance estimates of the
#' local component for the \eqn{j}-th individual in block \eqn{i} at time \eqn{t}.
#'
#' @param object An S3 object of class 'multi_result' created by multilevel().
#' @param i An integer indicating the \eqn{i}-th block.
#' @param j An integer indicating the \eqn{j}-th individual in the \eqn{i}-th block.
#' @param t An integer indicating the time.
#'
#' @return The variance of the local component.
#'
#' @export
#'
#' @examples
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' vcov_local_comp_ijt <- vcov_local_comp(est_multi, i = 1, j = 1, t = 1)
vcov_local_comp <- function(object, i, j, t) {
  G <- object$G
  F <- object$F
  Gamma <- object$Gamma
  Lambda <- object$Lambda
  R <- object$R
  N <- object$N
  Ni <- object$Ni
  T <- object$T
  r0 <- object$r0
  ri <- object$ri
  Resid <- object$Resid
  Y_block <- object$Y_list
  block <- i
  ind <- j

  Yi_hat <- Y_block[[block]] - G %*% t(Gamma[[block]])
  eig <- eigen(Yi_hat %*% t(Yi_hat) / (Ni[block] * T))
  Upsi_hat <- matrix(0, ri[block], ri[block])
  diag(Upsi_hat) <- eig$values[1:ri[block]]

  D3ij <- vcov_local_loading(object, i = block, j = ind)
  Var_Fi <- Upsi_hat %*% vcov_local_factor(object, i = block, t = t) %*% Upsi_hat

  Lambda_inv <- solve(t(Lambda[[block]]) %*% Lambda[[block]] / Ni[block])
  lambda_ij <- matrix(Lambda[[block]][ind, ], ncol = 1)
  Fit <- matrix(F[[block]][t, ], ncol = 1)

  Var_lc_ijt <- (t(lambda_ij) %*% Lambda_inv %*% Var_Fi %*% Lambda_inv %*% lambda_ij / Ni[block] + t(Fit) %*% D3ij %*% Fit / T) * min(Ni[block], T)

  return(c(Var_lc_ijt))
}

#' POET estimation of the covariance for the error terms
#'
#' @description
#' This function generates POET covariance estimates for the error terms.
#'
#' @param object An S3 object of class 'multi_result' created by [multilevel()].
#' @param C A positive constant in the adaptive threshold.
#'
#' @return An \eqn{N \times N} covariance matrix.
#' @import Matrix
#' @export
#'
#' @examples
#'
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' Sigma_e_POET <- vcov_POET(est_multi)
vcov_POET <- function(object, C = 1.5) {

  G <- object$G
  F <- object$F
  Gamma <- object$Gamma
  Lambda <- object$Lambda
  R <- object$R
  N <- object$N
  Ni <- object$Ni
  T <- object$T
  r0 <- object$r0
  ri <- object$ri
  Resid <- object$Resid

  e_mat <- matrix(unlist(Resid), nrow = T)
  Sigma_hat <- t(e_mat) %*% e_mat / T
  vSigma_hat <- t(e_mat^2) %*% (e_mat^2) / T - Sigma_hat^2
  Sigma_hat_thres <- sign(Sigma_hat) * pmax(abs(Sigma_hat) - C * sqrt(log(min(Ni)) / T) * sqrt(vSigma_hat), 0)
  rm(vSigma_hat)
  diag(Sigma_hat_thres) <- diag(Sigma_hat)

  return(Sigma_hat_thres)
}
