
#' Get bootstrap confidence intervals for the global factors
#'
#' @description
#' This function employs a bootstrap procedure to obtain confidence intervals
#' for the global factors at time \eqn{t}.
#'
#' @param object An S3 object of class 'multi_result' created by multilevel().
#' @param t An integer specifying the time point at which the CI is constructed.
#' @param BB An integer indicating the number of bootstrap repetition. 599 by default.
#' @param alpha The significance level, a single numeric between 0 and 1. 0.05 by default.
#'
#' @return A matrix containing the upper and lower band.
#' @export
#'
#' @examples
#'
#' \donttest{
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' bs_global_mid <- BS_global_factor(est_multi, t = est_multi$T / 2)
#' }
BS_global_factor <- function(object, t, BB = 599, alpha = 0.05) {

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

  bs_global <- sapply(c(1:BB), function(bb) {
    Resid_star <- lapply(c(1:R), function(i) {
      Resid[[i]] * matrix(stats::rnorm(T * Ni[i], 0, 1), T, Ni[i])
    })

    Y_block_star <- lapply(c(1:R), function(i) {
      if (ri[i] > 0) {
        G %*% t(Gamma[[i]]) + F[[i]] %*% t(Lambda[[i]]) + Resid_star[[i]]
      } else {
        G %*% t(Gamma[[i]]) + Resid_star[[i]]
      }
    })
    est_star <- multilevel(Y_block_star, r0 = r0, ri = ri)
    G_star <- est_star$G
    F_star <- est_star$F

    HG_star <- (1 / T) * solve(t(G) %*% G / T) %*% t(G) %*% G_star

    B_star_hat <- matrix(0, r0, r0)

    for (i in 1:R) {
      if (ri[i] > 0) {
        Ii <- rbind(diag(1, r0, r0), matrix(0, ri[i], r0))
        Thetai <- cbind(Gamma[[i]], Lambda[[i]])
        Ki_hat <- cbind(G, F[[i]])
        Pc <- PC(Y_block_star[[i]], r0 + ri[i])
        Ki_star_hat <- Pc$factor
        eig <- eigen(Y_block_star[[i]] %*% t(Y_block_star[[i]]) / (Ni[i] * T))
        Vi_star <- diag(eig$values[1:(r0 + ri[i])], r0 + ri[i], r0 + ri[i])
      } else {
        Ii <- rbind(diag(1, r0, r0))
        Thetai <- Gamma[[i]]
        Ki_hat <- G
        Pc <- PC(Y_block_star[[i]], r0)
        Ki_star_hat <- Pc$factor
        eig <- eigen(Y_block_star[[i]] %*% t(Y_block_star[[i]]) / (Ni[i] * T))
        Vi_star <- diag(eig$values[1:r0], r0, r0)
      }

      Hi_star <- (t(Thetai) %*% Thetai / Ni[i]) %*% (t(Ki_hat) %*% Ki_star_hat / T) %*% solve(Vi_star)

      U_star <- diag(diag(sign(t(G) %*% G_star / T)), r0, r0)
      B_star_hat <- B_star_hat + (1 / R) * (1 / Ni[i]) * t(Ii) %*% t(solve(Hi_star)) %*% solve(Vi_star) %*% (t(Ki_star_hat) %*% Ki_hat / T) %*%
        (t(Thetai) %*% t(Resid_star[[i]]) / (sqrt(Ni[i] * T))) %*% (G / sqrt(T)) %*% U_star
    }

    G_bs_t <- sqrt(N) * (solve(t(HG_star + B_star_hat)) %*% matrix(G_star[t, ], ncol = 1) - matrix(G[t, ], ncol = 1))

    return(G_bs_t)
  })

  bs_global <- matrix(t(bs_global), ncol = r0)

  upper <- G[t, ] + apply(bs_global, 2, stats::quantile, (1 - alpha) + alpha / 2) / sqrt(N)
  lower <- G[t, ] + apply(bs_global, 2, stats::quantile, alpha / 2) / sqrt(N)

  ret_mat <- matrix(c(upper, lower), ncol = 2)
  colnames(ret_mat) <- c("upper", "lower")

  return(ret_mat)
}

#' Get a bootstrap confidence interval for the global factor loadings
#'
#' @description This function employs a bootstrap procedure to obtain confidence intervals
#' for the global factor loadings for the \eqn{j}-th individual in block \eqn{i}. See Lin
#' and Shin (2023) for details.
#'
#' @param object An S3 object of class 'multi_result' created by [multilevel()].
#' @param i An integer indicating the \eqn{i}-th block.
#' @param j An integer indicating the \eqn{j}-th individual in the \eqn{i}-th block.
#' @param BB An integer indicating the number of bootstrap repetition. 599 by default.
#' @param alpha The significance level, a single numeric between 0 and 1. 0.05 by default.
#'
#' @return A matrix containing the upper and lower band.
#'
#' @references Lin, R. and Shin, Y., 2022. Generalised Canonical Correlation Estimation
#' of the Multilevel Factor Model. Available at SSRN 4295429.
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' bs_gamma_11 <- BS_global_loading(est_multi, i = 1, j = 1)
#' }
BS_global_loading <- function(object, i, j, BB = 599, alpha = 0.05) {
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
  Gamma_ij <- Gamma[[block]][ind, ]

  bs_global_loading <- sapply(c(1:BB), function(bb){
    Resid_star <- lapply(c(1:R), function(i) {
      Resid[[i]] * matrix(stats::rnorm(T * Ni[i], 0, 1), T, Ni[i])
    })

    Y_block_star <- lapply(c(1:R), function(i) {
      if (ri[i] > 0) {
        G %*% t(Gamma[[i]]) + F[[i]] %*% t(Lambda[[i]]) + Resid_star[[i]]
      } else {
        G %*% t(Gamma[[i]]) + Resid_star[[i]]
      }
    })

    FF <- lapply(c(1:R), function(i) {
      if (ri[i] > 0 & i == block) {
        F_mat <- matrix(0, T, ri[i])

        for (k in 1:ri[i]) {
          F_mat[, k] <- dwBS(matrix(F[[i]][, k], ncol = 1))
        }

        F_mat
      } else {
        F[[i]]
      }
    })

    Y_block_star <- lapply(c(1:R), function(i) {
      if (ri[i] > 0) {
        G %*% t(Gamma[[i]]) + FF[[i]] %*% t(Lambda[[i]]) + Resid_star[[i]]
      } else {
        G %*% t(Gamma[[i]]) + Resid_star[[i]]
      }
    })

    est_star <- multilevel(Y_block_star, r0 = r0, ri = ri)
    G_star <- est_star$G
    F_star <- est_star$F
    Gamma_star <- est_star$Gamma
    Lambda_star <- est_star$Lambda
    Gamma_star_ij <- Gamma_star[[block]][ind, ]

    HG_star <- (1 / T) * solve(t(G) %*% G / T) %*% t(G) %*% G_star
    B_star_hat <- matrix(0, r0, r0)

    for (i in 1:R) {
      if (ri[i] > 0) {
        Ii <- rbind(diag(1, r0, r0), matrix(0, ri[i], r0))
        Thetai <- cbind(Gamma[[i]], Lambda[[i]])
        Ki_hat <- cbind(G, FF[[i]])
        Pc <- PC(Y_block_star[[i]], r0 + ri[i])
        Ki_star_hat <- Pc$factor
        eig <- eigen(Y_block_star[[i]] %*% t(Y_block_star[[i]]) / (Ni[i] * T))
        Vi_star <- diag(eig$values[1:(r0 + ri[i])], r0 + ri[i], r0 + ri[i])
      } else {
        Ii <- diag(1, r0, r0)
        Thetai <- Gamma[[i]]
        Ki_hat <- G
        Pc <- PC(Y_block_star[[i]], r0)
        Ki_star_hat <- Pc$factor
        eig <- eigen(Y_block_star[[i]] %*% t(Y_block_star[[i]]) / (Ni[i] * T))
        Vi_star <- diag(eig$values[1:r0], r0, r0)
      }

      Hi_star <- (t(Thetai) %*% Thetai / Ni[i]) %*% (t(Ki_hat) %*% Ki_star_hat / T) %*% solve(Vi_star)

      U_star <- diag(diag(sign(t(G) %*% G_star / T)), r0, r0)
      B_star_hat <- B_star_hat + (1 / R) * (1 / Ni[i]) * t(Ii) %*% t(solve(Hi_star)) %*% solve(Vi_star) %*% (t(Ki_star_hat) %*% Ki_hat / T) %*%
        (t(Thetai) %*% t(Resid_star[[i]]) / (sqrt(Ni[i] * T))) %*% (G / sqrt(T)) %*% U_star
    }

    Gamma_bs_ij <- sqrt(T) * ((HG_star + B_star_hat) %*% matrix(Gamma_star_ij, ncol = 1)
      - matrix(Gamma_ij, ncol = 1))
    return(Gamma_bs_ij)
  })

  bs_global_loading <- matrix(t(bs_global_loading), ncol = r0)
  upper <- Gamma_ij + apply(bs_global_loading, 2, stats::quantile, 1 - alpha / 2) / sqrt(T)
  lower <- Gamma_ij + apply(bs_global_loading, 2, stats::quantile, alpha / 2) / sqrt(T)

  ret_mat <- matrix(c(upper, lower), ncol = 2)
  colnames(ret_mat) <- c("upper", "lower")

  return(ret_mat)
}


#' Get a bootstrap confidence interval for the global component
#'
#' @description This function employs a bootstrap procedure to obtain confidence intervals
#' for the global component for the \eqn{j}-th individual in block \eqn{i} at time \eqn{t}.
#' See Lin and Shin (2023) for details.
#'
#' @param object An S3 object of class 'multi_result' created by multilevel().
#' @param i An integer indicating the \eqn{i}-th block.
#' @param j An integer indicating the \eqn{j}-th individual in the \eqn{i}-th block.
#' @param t An integer specifying the time point at which the CI is constructed.
#' @param BB An integer indicating the number of bootstrap repetition. 599 by default.
#' @param alpha The significance level, a single numeric between 0 and 1. 0.05 by default.
#'
#' @return A matrix containing the upper and lower band.
#'
#' @references Lin, R. and Shin, Y., 2022. Generalised Canonical Correlation Estimation
#' of the Multilevel Factor Model. Available at SSRN 4295429.
#'
#' @export
#'
#' @examples
#' \donttest{
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' bs_gcomp_111 <- BS_global_comp(est_multi, i = 1, j = 1, t = 1)
#' }
BS_global_comp <- function(object, i, j, t, BB = 599, alpha = 0.05) {
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

  bs_global_comp <- sapply(c(1:BB), function(bb) {
    Resid_star <- lapply(1:R, function(i) {
      Resid[[i]] * matrix(stats::rnorm(T * Ni[i], 0, 1), T, Ni[i])
    })

    FF <- lapply(c(1:R), function(i) {
      if (ri[i] > 0 & i == block) {
        F_mat <- matrix(0, T, ri[i])

        for (k in 1:ri[i]) {
          F_mat[, k] <- dwBS(matrix(F[[i]][, k], ncol = 1))
        }

        F_mat
      } else {
        F[[i]]
      }
    })

    Y_block_star <- lapply(c(1:R), function(i) {
      if (ri[i] > 0) {
        G %*% t(Gamma[[i]]) + FF[[i]] %*% t(Lambda[[i]]) + Resid_star[[i]]
      } else {
        G %*% t(Gamma[[i]]) + Resid_star[[i]]
      }
    })
    est_star <- multilevel(Y_block_star, r0 = r0, ri = ri)
    G_star <- est_star$G
    G_t_star <- G_star[t, ]
    Gamma_star <- est_star$Gamma
    Gamma_ij_star <- Gamma_star[[block]][ind, ]
    Gcomp_ijt_star <- t(matrix(Gamma_ij_star, ncol = 1)) %*% matrix(G_t_star, ncol = 1)

    return(c(Gcomp_ijt_star))
  })

  upper <- stats::quantile(bs_global_comp, 1 - alpha / 2)
  lower <- stats::quantile(bs_global_comp, alpha/2)
  ret_mat <- matrix(c(upper, lower), ncol = 2)
  colnames(ret_mat) <- c("upper", "lower")

  return(ret_mat)
}

#' Get a bootstrap confidence interval for the local factors
#'
#' @description This function employs a bootstrap procedure to obtain confidence intervals
#' for the local factors in block \eqn{i} at time \eqn{t}. See Lin and Shin (2023) for details.
#'
#' @param object An S3 object of class 'multi_result' created by multilevel().
#' @param i An integer indicating the \eqn{i}-th block.
#' @param t An integer specifying the time point at which the CI is constructed.
#' @param BB An integer indicating the number of bootstrap repetition. 599 by default.
#' @param alpha The significance level, a single numeric between 0 and 1. 0.05 by default.
#'
#' @return A matrix containing the upper and lower band.
#'
#' @references Lin, R. and Shin, Y., 2022. Generalised Canonical Correlation Estimation
#' of the Multilevel Factor Model. Available at SSRN 4295429.
#'
#' @export
#'
#' @examples
#' \donttest{
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' bs_local_factor_11 <- BS_local_factor(est_multi, i = 1, t = 1)
#' }
BS_local_factor <- function(object, i, t, BB = 599, alpha = 0.05) {
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
  Fi_hat_t <- F[[block]][t, ]

  bs_local <- sapply(c(1:BB), function(bb) {
    G_new <- matrix(0, T, r0)

    for (k in 1:r0) {
      G_new[, k] <- dwBS(matrix(G[, k], ncol = 1))
    }

    Resid_star <- lapply(1:R, function(i) {
      Resid[[i]] * matrix(stats::rnorm(T * Ni[i], 0, 1), T, Ni[i])
    })
    Y_block_star <- lapply(c(1:R), function(i) {
      if (ri[i] > 0) {
        G_new %*% t(Gamma[[i]]) + F[[i]] %*% t(Lambda[[i]]) + Resid_star[[i]]
      } else {
        G_new %*% t(Gamma[[i]]) + Resid_star[[i]]
      }
    })
    est_star <- multilevel(Y_block_star, r0 = r0, ri = ri)
    G_star <- est_star$G
    F_star <- est_star$F
    Gamma_star <- est_star$Gamma
    Lambda_star <- est_star$Lambda

    Yi_star_hat <- Y_block_star[[block]] - G_star %*% t(Gamma_star[[block]])
    eig <- eigen((1 / (Ni[block] * T)) * Yi_star_hat %*% t(Yi_star_hat))
    Upsi_star_hat <- diag(eig$values[1:ri[block]], ri[block])

    F_bs_t <- matrix(F_star[[block]][t, ], ncol = 1)

    Hi_star <- (t(Lambda[[block]]) %*% Lambda[[block]] / Ni[block]) %*%
      (t(F[[block]]) %*% F_star[[block]] / T) %*% solve(Upsi_star_hat)


    Si_hat_star <- G_new %*% t(Gamma[[block]]) - G_star %*% t(Gamma_star[[block]])
    Bit_star <- solve(Upsi_star_hat) %*% (t(F_star[[block]]) %*% F[[block]] / T) %*%
      (t(Lambda[[block]]) %*% matrix(Si_hat_star[t, ], ncol = 1) / Ni[block])

    F_bs_t_rotated <- solve(t(Hi_star)) %*% F_bs_t
    diff <- F_bs_t_rotated - Fi_hat_t

    return(diff)
  })

  bs_local <- matrix(t(bs_local), ncol = ri[block])

  upper <- Fi_hat_t + apply(bs_local, 2, stats::quantile, (1 - alpha) + alpha / 2)
  lower <- Fi_hat_t + apply(bs_local, 2, stats::quantile, alpha / 2)
  ret_mat <- matrix(c(upper, lower), ncol = 2)
  colnames(ret_mat) <- c("upper", "lower")

  return(ret_mat)
}

#' Get a bootstrap confidence interval for the global component
#'
#' @description This function employs a bootstrap procedure to obtain confidence intervals
#' for the local component for the \eqn{j}-th individual in block \eqn{i} at time \eqn{t}.
#' See Lin and Shin (2023) for details.
#'
#' @param object An S3 object of class 'multi_result' created by multilevel().
#' @param i An integer indicating the \eqn{i}-th block.
#' @param j An integer indicating the \eqn{j}-th individual in the \eqn{i}-th block.
#' @param t An integer specifying the time point at which the CI is constructed.
#' @param BB An integer indicating the number of bootstrap repetition. 599 by default.
#' @param alpha The significance level, a single numeric between 0 and 1. 0.05 by default.
#'
#' @return A matrix containing the upper and lower band.
#'
#' @references Lin, R. and Shin, Y., 2022. Generalised Canonical Correlation Estimation
#' of the Multilevel Factor Model. Available at SSRN 4295429.
#'
#' @export
#'
#' @examples
#' \donttest{
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' bs_fcomp_111 <- BS_local_comp(est_multi, i = 1, j = 1, t = 1)
#' }
BS_local_comp <- function(object, i, j, t, BB = 599, alpha = 0.05) {
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

  bs_local_comp <- sapply(c(1:BB), function(bb) {
    G_new <- matrix(0, T, r0)
    for (k in 1:r0) {
      G_new[, k] <- dwBS(matrix(G[, k], ncol = 1))
    }


    Resid_star <- lapply(1:R, function(i) {
      Resid[[i]] * matrix(stats::rnorm(T * Ni[i], 0, 1), T, Ni[i])
    })
    Y_block_star <- lapply(c(1:R), function(i) {
      if (ri[i] > 0) {
        G_new %*% t(Gamma[[i]]) + F[[i]] %*% t(Lambda[[i]]) + Resid_star[[i]]
      } else {
        G_new %*% t(Gamma[[i]]) + Resid_star[[i]]
      }
    })
    est_star <- multilevel(Y_block_star, r0 = r0, ri = ri)
    G_star <- est_star$G
    F_star <- est_star$F
    Gamma_star <- est_star$Gamma
    Lambda_star <- est_star$Lambda
    F_t_star <- F_star[[block]][t, ]
    Lambda_ij_star <- Lambda_star[[block]][ind, ]

    Fcomp_bs_ijt <- t(matrix(Lambda_ij_star, ncol = 1)) %*%
      matrix(F_t_star, ncol = 1)

    return(Fcomp_bs_ijt)
  })

  bs_local_comp <- t(bs_local_comp)
  upper <- stats::quantile(bs_local_comp, 1 - alpha / 2)
  lower <- stats::quantile(bs_local_comp, alpha / 2)
  ret_mat <- matrix(c(upper, lower), ncol = 2)
  colnames(ret_mat) <- c("upper", "lower")

  return(ret_mat)
}

#' Get an asymptotic confidence interval for the local component
#'
#' @description This function computes the asymptotic confidence intervals
#' for the local loadings for the \eqn{j}-th individual in block \eqn{i}.
#' See Lin and Shin (2023) for details.
#'
#' @param object An S3 object of class 'multi_result' created by multilevel().
#' @param i An integer indicating the \eqn{i}-th block.
#' @param j An integer indicating the \eqn{j}-th individual in the \eqn{i}-th block.
#' @param alpha The significance level, a single numeric between 0 and 1. 0.05 by default.
#'
#' @return A matrix containing the upper and lower band.
#'
#' @references Lin, R. and Shin, Y., 2022. Generalised Canonical Correlation Estimation
#' of the Multilevel Factor Model. Available at SSRN 4295429.
#'
#' @export
#'
#' @examples
#' \donttest{
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' bs_local_loading_11 <- AsymCI_local_loading(est_multi, i = 1, j = 1)
#' }
AsymCI_local_loading <- function(object, i, j, alpha = 0.05) {

  F <- object$F
  Lambda <- object$Lambda
  T <- object$T
  Resid <- object$Resid
  Lambda_ij_hat <- Lambda[[i]][j, ]
  FEij <- F[[i]] * c(Resid[[i]][, j])
  fit <- stats::lm(FEij ~ 1)

  Var_FEij <- sandwich::NeweyWest(fit, sandwich = FALSE)

  Var_Lambdaij_hat <- Var_FEij

  upper <- Lambda_ij_hat + sqrt(diag(Var_Lambdaij_hat)) * stats::qnorm(1 - alpha / 2) / sqrt(T)
  lower <- Lambda_ij_hat - sqrt(diag(Var_Lambdaij_hat)) * stats::qnorm(1 - alpha / 2) / sqrt(T)

  ret_mat <- matrix(c(upper, lower), ncol = 2)
  colnames(ret_mat) <- c("upper", "lower")

  return(ret_mat)
}
