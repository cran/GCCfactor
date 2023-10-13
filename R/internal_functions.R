#' Trace of a matrix
#' @description
#' This function calculates the trace of a square matrix \eqn{tr(\boldsymbol{A}).}
#'
#' @param m A square matrix.
#'
#' @return The trace of the matrix.
#' @noRd
tr <- function(m) {
  total_sum <- 0
  if (is.matrix(m)) {
    row_count <- nrow(m)
    col_count <- ncol(m)
    if (row_count == col_count) {
      total_sum <- sum(diag(m))
      total_sum
    } else {
      message("Matrix is not square")
    }
  } else {
    message("Object is not a matrix")
  }
}

#' Get the matrix \eqn{\widehat{\boldsymbol{\Phi}}} which is used in GCC().
#'
#' @param K A list containing the estimated factors \eqn{\widehat{\boldsymbol{K}}_{i}} extracted from each block.
#' @param R An integer of the number of blocks.
#' @param d An array of length \eqn{R} containing the total number of factors in each block.
#' @param T = number of time series observations.
#'
#' @return The matrix \eqn{\widehat{\boldsymbol{\Phi}}}.
#' @noRd
get_Phi <- function(K, R, d, T) {
  Phi <- matrix(0, nrow = T * R * (R - 1) / 2, ncol = sum(d))
  rowstart <- 1
  rowend <- T
  csumd <- cumsum(d)
  for (m in 1:(R - 1)) {
    for (h in (m + 1):R) {
      if (m == 1) {
        mstart <- 1
      } else {
        mstart <- csumd[m - 1] + 1
      }
      mend <- csumd[m]

      hstart <- csumd[h - 1] + 1
      hend <- csumd[h]

      Phi[rowstart:rowend, mstart:mend] <- K[[m]]
      Phi[rowstart:rowend, hstart:hend] <- -K[[h]]

      rowstart <- rowstart + T
      rowend <- rowend + T
    }
  }
  Phi
}

#' Supply numbers of factors in each block.
#'
#' @param r_max An integer indicating the maximum number of factors allowed in each block.
#' @param r0 An integer of the number of global factors.
#' @param ri An array of length \eqn{R} containing the number of local factors in each block.
#' @param R An integer of the number of blocks.
#'
#' @return An array of length \eqn{R} containing the total number of factors in each block.
#' @noRd
get_K_dim <- function(r_max, r0, ri, R) {
  r0_empty <- is.null(r0)
  ri_empty <- is.null(ri)
  if (r0_empty | ri_empty) {
    if((r0_empty & !ri_empty) | (!r0_empty & ri_empty)){
      message(paste0("r0 or ri is NULL. The one that is NULL will be estimated."))
    }
    d <- rep(r_max, R)
  }else if (!r0_empty & !ri_empty) {
    if (length(ri) != R) {
      stop(paste0("Length of ri should be equal to the number of blocks."))
    }
    d <- r0 + ri
  }else {
    stop("Numbers of factors are not correctly specified. Try r0 = NULL and ri = NULL.")
  }
  return(d)
}

