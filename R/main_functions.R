#' Generalised canonical correlation estimation for the global factors
#'
#' @description
#' This function is one of the main functions the package, employing the
#' generalized canonical correlation estimation for both the global factors
#' \eqn{\boldsymbol{G}} and, when not explicitly provided, for the number of
#' global factors \eqn{r_{0}}. Typically, this function is intended for internal
#' purposes. However, users one can opt for GCC() instead of multilevel(),
#' if the users only need to estimate the number of global factors.
#'
#'
#' @param data Either a data.frame or a list of data matrices of length \eqn{R}. See \strong{Details}.
#' @param standarise A logical indicating whether the data is standardised before estimation or not. See \strong{Details}.
#' @param r_max An integer indicating the maximum number of factors allowed. See \strong{Details}.
#' @param r0 An integer of the number of global factors. See \strong{Details}.
#' @param ri An array of length \eqn{R} containing the number of local factors in each block. See \strong{Details}.
#' @param depvar_header A character string specifying the header of the dependent variable. See \strong{Details}.
#' @param i_header A character string specifying the header of the block identifier. See \strong{Details}.
#' @param j_header A character string specifying the header of the individual identifier. See \strong{Details}.
#' @param t_header A character string specifying the header of the time identifier. See \strong{Details}.
#'
#' @details
#' The user-supplied data.frame should contain at least four columns, namely the
#' dependent variable (\eqn{y_{ijt}}), block identifier (\eqn{i}), individual
#' identifier (\eqn{j}), and time (\eqn{t}). The user needs to supply their corresponding
#' headers in the data.frame to the function using the parameters "depvar_header",
#' "i_header", "j_header", and "t_header", respectively. If the data is supplied
#' as a list, these arguments will not be used.
#'
#' If either r0 = NULL or ri = NULL, both of them will be estimated.
#' In such case, "r_max" must be supplied. If "r0" and "ri" are supplied then
#' "r_max" is not needed and will be ignored.
#'
#' If standarise = TRUE, each time series will be standardised so it has zero mean
#' and unit variance. It is recommended to standardise the data before estimation.
#'
#' See Lin and Shin (2023) for more details.
#'
#' @return A list containing the estimated number of global factors \eqn{\hat{r}_{0}},
#' the global factors \eqn{\widehat{\boldsymbol{G}}}, and the other elements that are
#' used in multilevel().
#'
#' @references Lin, R. and Shin, Y., 2022. Generalised Canonical Correlation Estimation
#' of the Multilevel Factor Model. Available at SSRN 4295429.
#'
#' @export
#'
#' @examples
#'
#' panel <- UKhouse # load the data
#' Y_list <- panel2list(panel, depvar_header = "dlPrice", i_header = "Region",
#'                                        j_header = "LPA_Type", t_header = "Date")
#' est_GCC <- GCC(Y_list, r_max = 10)
#' r0_hat <- est_GCC$r0 # number of global factors
#' G_hat <- est_GCC$G # global factors
GCC <- function(data, standarise = TRUE, r_max = 10, r0 = NULL, ri = NULL, depvar_header = NULL,
                i_header = NULL, j_header = NULL, t_header = NULL) {
  Y_list <- check_data(data,
    depvar_header = depvar_header, i_header = i_header,
    j_header = j_header, t_header = t_header
  )

  if (isTRUE(standarise)) {
    Y_list <- lapply(Y_list, scale)
  }

  R <- length(Y_list)
  Ni <- sapply(Y_list, ncol)
  T <- nrow(Y_list[[1]])

  d <- get_K_dim(r_max, r0, ri, R)
  K <- lapply(c(1:R), function(i) {
    est <- PC(Y_list[[i]], d[i])
    return(est$factor)
  })

  Phi <- get_Phi(K, R, d, T)
  s <- svd(Phi / sqrt(T))
  v <- s$v
  delta2 <- (s$d[length(s$d):1])^2
  mock <- mean(delta2) / (sqrt(min(T, Ni)))
  delta2 <- delta2[1:(r_max + 1)]
  delta2 <- c(mock, delta2)

  if (is.null(r0)) {
    ratios <- delta2[2:length(delta2)] / delta2[1:(length(delta2) - 1)]
    r0 <- which(ratios == max(ratios)) - 1
  }

  if (r0 > 0) {
    Q <- matrix(v[, (ncol(v) - r0 + 1):ncol(v)], nrow = nrow(v))
    csumd <- cumsum(d)
    Psi <- lapply(1:R, function(i) {
      if (i == 1) {
        rowstart <- 1
      } else {
        rowstart <- csumd[i - 1] + 1
      }
      rowend <- csumd[i]
      K[[i]] %*% Q[rowstart:rowend, ]
    })

    Psi <- matrix(unlist(Psi), nrow = T)
    s <- svd(Psi)
    G <- sqrt(T) * matrix(s$u[, 1:r0], nrow = T)
  } else if (r0 == 0) {
    G <- NA
  }
  return_val <- list(G, r0, d, R, Ni, T, delta2, Y_list)
  names(return_val) <- c("G", "r0", "d", "R", "Ni", "T", "delta2", "Y_list")
  return(return_val)
}


#' Full estimation of the multilevel factor model
#'
#' @description
#' This is one of the main functions of this package which performs full estimation
#' of the multilevel factor model.
#'
#'
#' @param data Either a data.frame or a list of data matrices of length \eqn{R}. See \strong{Details}.
#' @param ic A character string of selection criteria to use for estimation of the numbers of local factors. See \strong{Details}.
#' @param standarise A logical indicating whether the data is standardised before estimation or not. See \strong{Details}.
#' @param r_max An integer indicating the maximum number of factors allowed. See \strong{Details}.
#' @param r0 An integer of the number of global factors. See \strong{Details}.
#' @param ri An array of length \eqn{R} containing the number of local factors in each block. See \strong{Details}.
#' @param depvar_header A character string specifying the header of the dependent variable. See \strong{Details}.
#' @param i_header A character string specifying the header of the block identifier. See \strong{Details}.
#' @param j_header A character string specifying the header of the individual identifier. See \strong{Details}.
#' @param t_header A character string specifying the header of the time identifier. See \strong{Details}.
#'
#' @details
#' The user-supplied data.frame should contain at least four columns, namely the
#' dependent variable (\eqn{y_{ijt}}), block identifier (\eqn{i}), individual
#' identifier (\eqn{j}), and time (\eqn{t}). The user needs to supply their corresponding
#' headers in the data.frame to the function using the parameters "depvar_header",
#' "i_header", "j_header", and "t_header", respectively. If the data is supplied
#' as a list, these arguments will not be used.
#'
#' If either r0 = NULL or ri = NULL, then both of them will be estimated.
#' In such case, "r_max" must be supplied. If "r0" and "ri" are supplied then
#' "r_max" is not needed and will be ignored.
#'
#' If standarise = TRUE, each time series will be standardised so it has zero mean
#' and unit variance. It is recommended to standardise the data before estimation.
#'
#' See Lin and Shin (2023) for more details.
#'
#'
#' @return The return value is an S3 object of class "multi_result".
#' It contains a list of the following items:
#' \itemize{
#'
#' \item G = A matrix of the estimated global factors.
#' \item Gamma = A list of length \eqn{R} containing matrices of the estimated global loading matrices for each block.
#' \item F = A list of length \eqn{R} containing matrices of the estimated local factors for each block.
#' \item Lambda = A list of length \eqn{R} containing matrices of the estimated global loading matrices for each block.
#' \item N = The total number of cross-sections in the panel.
#' \item Ni = An array of length \eqn{R} containing the number of cross-sections in each block.
#' \item r0 = The number of global factors. Unchanged if pre-specified.
#' \item ri = An array of length \eqn{R} containing the number of local factors for each block. Unchanged if pre-specified.
#' \item d = An array of length \eqn{R} containing the maximum total number of factors allowed for each block.
#' The elements are identically equal to r_max if either r0 or ri is supplied as NULL.
#' \item Resid = A list of length \eqn{R} containing the residual matrices for each block.
#' \item delta2 = An array of the mock and the \eqn{r_{\max} + 1} largest squared singular values.
#' \item ic = Selection criteria used for estimating the numbers of local factors.
#' \item block_names = A array of block names.
#'
#' }
#'
#' @references Lin, R. and Shin, Y., 2022. Generalised Canonical Correlation Estimation
#' of the Multilevel Factor Model. Available at SSRN 4295429.
#'
#' @export
#'
#' @examples
#'
#' panel <- UKhouse # load the data
#'
#' # use data.frame
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' # or one can use a list of data matrices
#' Y_list <- panel2list(panel, depvar_header = "dlPrice", i_header = "Region",
#'                                        j_header = "LPA_Type", t_header = "Date")
#' est_multi <- multilevel(Y_list, ic = "BIC3", standarise = TRUE, r_max = 5)
multilevel <- function(data, ic = "BIC3", standarise = TRUE, r_max = 10, r0 = NULL,
                       ri = NULL, depvar_header = NULL, i_header = NULL, j_header = NULL, t_header = NULL) {
  GCC_est <- GCC(data, standarise, r_max, r0, ri, depvar_header, i_header, j_header, t_header)
  G <- GCC_est$G
  r0 <- GCC_est$r0
  d <- GCC_est$d
  R <- GCC_est$R
  Ni <- GCC_est$Ni
  delta2 <- GCC_est$delta2
  Y_list<- GCC_est$Y_list
  T <- GCC_est$T
  i_names <- names(Y_list)
  N <- sum(Ni)

  ri_max <- d - r0
  ri_empty <- is.null(ri)
  if (ri_empty) {
    ri <- c()
  }

  Gamma <- list()
  F <- list()
  Lambda <- list()
  Y_list_hat <- list()
  if (r0 > 0) {
    for (i in 1:R) {
      Gamma[[i]] <- (1 / T) * t(Y_list[[i]]) %*% G
      Y_list_hat[[i]] <- Y_list[[i]] - G %*% t(Gamma[[i]])
      if (ri_empty) {
        ri[i] <- infocrit(Y_list_hat[[i]], ic, ri_max[i])
      }

      if (ri[i] > 0) {
        est <- PC(Y_list_hat[[i]], ri[i])
        F[[i]] <- est$factor
        Lambda[[i]] <- est$loading
      } else {
        F[[i]] <- NA
        Lambda[[i]] <- NA
      }
    }
  } else {
    if (ri_empty) {
      ri <- sapply(Y_list, function(x) {
        infocrit(x, ic, ri_max[i])
      })
    }
    G <- NA
    Gamma <- NA

    for (i in 1:R) {
      if (ri[i] > 0) {
        est <- PC(Y_list[[i]], ri[i])
        F[[i]] <- est$factor
        Lambda[[i]] <- est$loading
      } else {
        F[[i]] <- NA
        Lambda[[i]] <- NA
      }
    }
  }

  Resid <- lapply(c(1:R), function(i) {
    if (r0 > 0 & ri[i] > 0) {
      resid <- Y_list[[i]] - G %*% t(Gamma[[i]]) - F[[i]] %*% t(Lambda[[i]])
      return(resid)
    } else if (r0 > 0 & ri[i] == 0) {
      resid <- Y_list[[i]] - G %*% t(Gamma[[i]])
      return(resid)
    } else if (r0 == 0 & ri[i] > 0) {
      resid <- Y_list[[i]] - F[[i]] %*% t(Lambda[[i]])
      return(resid)
    } else {
      resid <- Y_list[[i]]
      return(resid)
    }
  })

  return_val <- list(G, F, Gamma, Lambda, R, N, Ni, T, r0, ri, d, Resid, delta2, ic, i_names, Y_list)
  names(return_val) <- c("G", "F", "Gamma", "Lambda", "R", "N", "Ni", "T", "r0",
                         "ri", "d", "Resid", "delta2", "ic", "block_names", "Y_list")
  class(return_val) <- "multi_result"
  return(return_val)
}



