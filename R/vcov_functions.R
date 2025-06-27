
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
  for (i in 1:R) {
    Ei <- Resid[[i]]
    Ii <- rbind(diag(1, r0, r0), matrix(0, ri[i], r0))
    if(ri[i] > 0){
      Thetai <- cbind(Gamma[[i]], Lambda[[i]])
    }else{
      Thetai <- Gamma[[i]]
    }

    Diit <- matrix(0, r0 + ri[i], r0 + ri[i])
    for (j in 1:Ni[i]) {
      Diit <- Diit + (1 / Ni[i]) * ((Ei[t, j])^2) * matrix(Thetai[j, ], ncol = 1) %*% t(matrix(Thetai[j, ], ncol = 1))
    }

    Var_G <- Var_G + (R^(-2)) * (N / Ni[i]) * t(Ii) %*% solve(t(Thetai) %*% Thetai / Ni[i]) %*%
      Diit %*% solve(t(Thetai) %*% Thetai / Ni[i]) %*% Ii
  }

  return(Var_G)
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
#'
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
  Fi_hat_t <- matrix(F[[block]][t, ], ncol = 1)

  Ii <- rbind(diag(1, r0, r0), matrix(0, ri[block], r0))
  IIi <- rbind(matrix(0, r0, ri[block]), diag(1, ri[block], ri[block]))
  Yi_hat <- Y_block[[block]] - G %*% t(Gamma[[block]])
  eig <- eigen(Yi_hat %*% t(Yi_hat) / (Ni[block] * T))
  Upsi_hat <- diag(matrix(c(eig$values[1:ri[block]]), ncol = 1))

  Var_Fi <- 0
  for (m in 1:R) {
    if(ri[m] > 0){
      Thetam <- cbind(Gamma[[m]], Lambda[[m]])
    }else{
      Thetam <- Gamma[[m]]
    }

    Dmmt <- 0
    Ci <- 0
    for (j in 1:Ni[m]) {
      Dmmt <- Dmmt + (1 / Ni[m]) * ((Resid[[m]][t, j])^2) * matrix(Thetam[j, ], ncol = 1) %*% t(matrix(Thetam[j, ], ncol = 1))
    }
    if (m == block) {
      Ci <- t(IIi) - sqrt(Ni[block] / N) * sqrt(N / Ni[m]) * (t(Lambda[[block]]) %*% Gamma[[block]] / Ni[block]) %*% t(Ii) %*% solve(t(Thetam) %*% Thetam / Ni[m]) / R
      Var_Fi <- Var_Fi + Ci %*% Dmmt %*% t(Ci)
    } else {
      if(ri[m] > 0){
        Ci <- -sqrt(Ni[block] / N) * sqrt(N / Ni[m]) * (t(Lambda[[block]]) %*% Gamma[[block]] / Ni[block]) %*% t(Ii) %*% solve(t(Thetam) %*% Thetam / Ni[m]) / R
        Var_Fi <- Var_Fi + Ci %*% Dmmt %*% t(Ci)
      }else{
        Var_Fi <- Var_Fi
      }
    }
  }
  Var_Fi <- solve(Upsi_hat) %*% Var_Fi %*% solve(Upsi_hat)

  return(Var_Fi)
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
#'
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

  D3ij <- vcov_local_loading(object, i = block, j = ind)

  Ii <- rbind(diag(1, r0, r0), matrix(0, ri[block], r0))
  IIi <- rbind(matrix(0, r0, ri[block]), diag(1, ri[block], ri[block]))
  Yi_hat <- Y_block[[block]] - G %*% t(Gamma[[block]])

  Var_Fi <- 0
  for (m in 1:R) {
    if(ri[m] > 0){
      Thetam <- cbind(Gamma[[m]], Lambda[[m]])
    }else{
      Thetam <- Gamma[[m]]
    }

    Dmmt <- 0
    Ci <- 0
    for (j in 1:Ni[m]) {
      Dmmt <- Dmmt + (1 / Ni[m]) * ((Resid[[m]][t, j])^2) * matrix(Thetam[j, ], ncol = 1) %*% t(matrix(Thetam[j, ], ncol = 1))
    }
    if (m == block) {
      Ci <- t(IIi) - sqrt(Ni[block] / N) * sqrt(N / Ni[m]) * (t(Lambda[[block]]) %*% Gamma[[block]] / Ni[block]) %*% t(Ii) %*% solve(t(Thetam) %*% Thetam / Ni[m]) / R
      Var_Fi <- Var_Fi + Ci %*% Dmmt %*% t(Ci)
    } else {
      if(ri[m] > 0){
        Ci <- -sqrt(Ni[block] / N) * sqrt(N / Ni[m]) * (t(Lambda[[block]]) %*% Gamma[[block]] / Ni[block]) %*% t(Ii) %*% solve(t(Thetam) %*% Thetam / Ni[m]) / R
        Var_Fi <- Var_Fi + Ci %*% Dmmt %*% t(Ci)
      }else{
        Var_Fi <- Var_Fi
      }
    }
  }

  Lambda_inv <- solve(t(Lambda[[block]]) %*% Lambda[[block]] / Ni[block])
  lambda_ij <- matrix(Lambda[[block]][ind, ], ncol = 1)
  Fit <- matrix(F[[block]][t, ], ncol = 1)

  Var_lc_ijt <- (t(lambda_ij) %*% Lambda_inv %*% Var_Fi %*% Lambda_inv %*% lambda_ij / Ni[block] + t(Fit) %*% D3ij %*% Fit / T) * min(Ni[block], T)

  return(c(Var_lc_ijt))
}
