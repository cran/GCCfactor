#' Principal component (PC) estimation of the approximate factor model
#'
#' @description
#' Perform PC estimation of the (2D) approximate factor model:
#' \deqn{y_{it}=\boldsymbol{\lambda}_{i}^{\prime}\boldsymbol{F}_{t}+e_{it},}
#' or in matrix notation:
#' \deqn{\boldsymbol{Y}=\boldsymbol{F}\boldsymbol{\Lambda}^{\prime}+\boldsymbol{e}.}
#' The factors \eqn{\boldsymbol{F}} is estimated as \eqn{\sqrt{T}} times the \eqn{r} eigenvectors of
#' the matrix \eqn{\boldsymbol{Y}\boldsymbol{Y}^{\prime}} corresponding to the \eqn{r}
#' largest eigenvalues in descending order, and the loading matrix is estimated by
#' \eqn{\boldsymbol{\Lambda}=T^{-1}\boldsymbol{Y}^{\prime}\boldsymbol{F}}.
#' See e.g. Bai and Ng (2002).
#'
#' @param Y A \eqn{T \times N} data matrix. T = number of time series observations,
#' N = cross-sectional dimension.
#' @param r = the number of factors.
#'
#' @return A list containing the factors and factor loadings:
#' \itemize{
#'
#'    \item factor = a \eqn{T \times r} matrix of the estimated factors.
#'
#'    \item loading = a \eqn{N \times r} matrix of the estimated factor loadings.
#'
#' }
#'
#' @references Bai, J. and Ng, S., 2002. Determining the number of factors in approximate factor models.
#' Econometrica, 70(1), pp.191-221.
#'
#' @export
#'
#' @examples
#'
#' # simulate data
#'
#' T <- 100
#' N <- 50
#' r <- 2
#' F <- matrix(stats::rnorm(T * r, 0, 1), nrow = T)
#' Lambda <- matrix(stats::rnorm(N * r, 0, 1), nrow = N)
#' err <- matrix(stats::rnorm(T * N, 0, 1), nrow = T)
#' Y <- F %*% t(Lambda) + err
#'
#' # estimation
#'
#' est_PC <- PC(Y, r)
PC <- function(Y, r) {
  T <- nrow(Y)
  eig <- eigen(Y %*% t(Y))
  F <- matrix(sqrt(T) * eig$vectors[, 1:r], nrow = T, ncol = r)
  Lambda <- (1 / T) * t(Y) %*% F
  return_val <- list(F, Lambda)
  names(return_val) <- c("factor", "loading")
  return(return_val)
}

#' Selection criteria for the approximate factor model
#' @description
#' This function performs model selection for the (2D) approximate factor model
#' and returns the estimated number of factors.
#'
#' @param Y A \eqn{T \times N} data matrix. T = number of time series observations,
#' N = cross-sectional dimension.
#' @param method A character string indicating which criteria to use.
#'
#' @param r_max An integer indicating the maximum number of factors allowed. 10 by default.
#'
#' @details "method" can be one of the following: "ICp2" and "BIC3" by Bai and Ng (2002),
#' "ER" by Ahn and Horenstein (2013), "ED" by Onatski (2010).
#'
#' @return The estimated number of factors.
#'
#' @references
#' Bai, J. and Ng, S., 2002. Determining the number of factors in approximate factor models.
#' Econometrica, 70(1), pp.191-221.
#'
#' Ahn, S.C. and Horenstein, A.R., 2013. Eigenvalue ratio test for the number of factors.
#' Econometrica, 81(3), pp.1203-1227.
#'
#' Onatski, A., 2010. Determining the number of factors from empirical distribution of eigenvalues.
#' The Review of Economics and Statistics, 92(4), pp.1004-1016.
#'
#' @export
#'
#' @examples
#' # simulate data
#'
#' T <- 100
#' N <- 50
#' r <- 2
#' F <- matrix(stats::rnorm(T * r, 0, 1), nrow = T)
#' Lambda <- matrix(stats::rnorm(N * r, 0, 1), nrow = N)
#' err <- matrix(stats::rnorm(T * N, 0, 1), nrow = T)
#' Y <- F %*% t(Lambda) + err
#'
#' # estimation
#'
#' r_hat <- infocrit(Y, "BIC3", r_max = 10)
infocrit <- function(Y, method, r_max = 10) {
  N <- ncol(Y)
  T <- nrow(Y)

  if (method == "ICp2" | method == "BIC3") {
    # ICp2 and BIC3

    # demean the data
    Y <- scale(Y, scale = FALSE)
    est <- PC(Y, r_max)
    F_max <- est$factor
    L_max <- est$loading
    res_max <- Y - F_max %*% t(L_max)
    RSS_max <- (1 / (N * T)) * tr(res_max %*% t(res_max))
    ICp2 <- matrix(0, nrow = r_max + 1, ncol = 1)
    BIC3 <- matrix(0, nrow = r_max + 1, ncol = 1)
    for (r in 0:r_max)
    {
      if (r != 0) {
        est <- PC(Y, r)
        F <- est$factor
        L <- est$loading
        res <- Y - F %*% t(L)
        RSS <- (1 / (N * T)) * tr(res %*% t(res))
      } else {
        RSS <- (1 / (N * T)) * tr(Y %*% t(Y))
      }

      ICp2[r + 1, 1] <- log(RSS) + r * ((N + T) / (N * T)) * log(min(N, T))
      BIC3[r + 1, 1] <- RSS + r * RSS_max * ((N + T - r) * log(N * T) / (N * T))
    }
    if (method == "ICp2") {
      r_hat <- which(ICp2 == min(ICp2)) - 1
    } else if (method == "BIC3") {
      r_hat <- which(BIC3 == min(BIC3)) - 1
    }
  } else if (method == "ER") {
    # ER
    # Double demean
    X <- Y
    # Demean by CS
    ts_mean <- colMeans(Y)
    # Demean by TS
    cs_mean <- rowMeans(Y)
    # Add overall mean
    g_mean <- mean(Y)
    X <- X - matrix(rep(ts_mean, each = T), T, N) - matrix(rep(cs_mean, N), T, N) + g_mean

    D <- X %*% t(X) / (N * T)
    eig <- eigen(D)
    eval <- eig$values[1:(r_max + 1)]
    m <- min(N, T)
    mock <- sum(eig$values) / log(m)
    eval <- c(mock, eval)
    ER <- matrix(0, nrow = r_max + 1, ncol = 1)
    for (i in 1:(r_max + 1))
    {
      ER[i, 1] <- eval[i] / eval[i + 1]
    }
    r_hat <- which(ER == max(ER)) - 1
  } else if (method == "ED") {
    # ED
    # demean
    Y <- scale(Y, scale = FALSE)
    D <- stats::cov(Y)
    eig <- eigen(D)
    eval <- eig$values
    j <- r_max + 1
    for (count in 1:5) {
      lamb <- eval[j:(j + 4)]
      const <- c((j - 1)^(2 / 3), (j)^(2 / 3), (j + 1)^(2 / 3), (j + 2)^(2 / 3), (j + 3)^(2 / 3))
      reg <- stats::lm(lamb ~ const)
      b_hat <- as.numeric(reg$coefficients[2])
      delta <- 2 * abs(b_hat)
      ED <- eval[1:r_max] - eval[2:(r_max + 1)]
      if (any(ED > delta) == FALSE) {
        r_hat <- 0
        j <- r_hat + 1
      } else {
        r_hat <- max(which(ED > delta))
        j <- r_hat + 1
      }
    }
  }
  return(r_hat)
}
