#' data.frame to list of data matrices
#'
#' @description
#' This function converts the data.frame to a list of data matrices and finds
#' the dimensions of the multilevel panel.
#'
#' @param panel The user-supplied data frame for the multilevel panel data. See \strong{Details}.
#' @param depvar_header A character string specifying the header of the dependent variable. See \strong{Details}.
#' @param i_header A character string specifying the header of the block identifier. See \strong{Details}.
#' @param j_header A character string specifying the header of the individual identifier. See \strong{Details}.
#' @param t_header A character string specifying the header of the time identifier. See \strong{Details}.
#'
#' @details
#' See the details of GCC().
#'
#'
#' @return A list containing the data matrices of the \eqn{R} blocks. Each of them
#' has dimension \eqn{T\times N_{i}}.
#' @export
#'
#' @examples
#'
#' panel <- UKhouse # load the data
#'
#' # panel$Region identifies different blocks i=1,...,R.
#' # panel$LPA_Type identifies different individuals j=1,...,N_i.
#'
#' Y_list<- panel2list(panel, depvar_header = "dlPrice", i_header = "Region",
#'                                        j_header = "LPA_Type", t_header = "Date")
#'
panel2list <- function(panel, depvar_header = NULL, i_header = NULL,
                       j_header = NULL, t_header = NULL) {
  i_names <- unique(unlist(panel[, i_header]))
  R <- length(i_names)
  List <- lapply(c(1:R), function(i) {
    select <- panel[, i_header] == i_names[i]
    subpanel <- panel[select ,c(depvar_header, j_header, t_header)]
    subpanel_wide <- reshape2::dcast(subpanel,
                           formula = paste(t_header, "~", j_header),
                           value.var = depvar_header)
    subpanel_mat <- data.matrix(subpanel_wide)
    subpanel_mat[, 2:ncol(subpanel_mat)]
  })
  names(List) <- i_names
  return(List)
}

#' Check validity of the data and headers
#'
#' @description
#' This is an internal function which checks the validity of the data and
#' provide a list of matrices of length \eqn{R} for estimation.
#'
#' @param data Either a data.frame or a list of data matrices of length \eqn{R}. See \strong{Details}.
#' @param depvar_header A character string specifying the header of the dependent variable. See \strong{Details}.
#' @param i_header A character string specifying the header of the block identifier. See \strong{Details}.
#' @param j_header A character string specifying the header of the individual identifier. See \strong{Details}.
#' @param t_header A character string specifying the header of the time identifier. See \strong{Details}.
#'
#' @details
#' See \strong{Details} of [GCC()].
#'
#'
#' @return A list of data matrices of length \eqn{R}.
#' @export
#'
#' @examples
#' panel <- UKhouse # load the data
#' Y_list <- check_data(panel,
#'   depvar_header = "dlPrice", i_header = "Region",
#'   j_header = "LPA_Type", t_header = "Date"
#' )
check_data <- function(data, depvar_header = NULL, i_header = NULL,
                       j_header = NULL, t_header = NULL) {
  headers <- as.character(c(depvar_header, i_header, j_header, t_header))
  if (is.data.frame(data)) {
    if (length(headers) != 4) {
      stop("The header names are not correctly specified.")
    }
    Y_list <- panel2list(data, depvar_header, i_header,
                          j_header, t_header)
  } else if (is.list(data)) {
    if (length(headers) != 0) {
      message("A list of data is supplied. The header names are not used.")
    }
    if (length(names(data)) != length(data)) {
      names(data) <- paste0("Block ", c(1:length(data)))
    }
    Y_list <- data
  }
  nrows <- sapply(Y_list, nrow)
  if (!length(unique(nrows)) == 1) {
    stop("The panel data is not balanced.")
  }
  if (any(sapply(Y_list, function(x) {
    any(is.na(x))
  }))) {
    stop("Data contains NAs.")
  }
  return(Y_list)
}

#' Print the relative importance ratios
#'
#' @param object An S3 object of class 'multi_result' created by multilevel().
#' @param ... Additional arguments.
#'
#' @return A matrix containing the summary of the model.
#' @usage \method{summary}{multi_result}(object, ...)
#'
#' @examples
#'
#' panel <- UKhouse # load the data
#' est_multi <- multilevel(panel, ic = "BIC3", standarise = TRUE, r_max = 5,
#'                            depvar_header = "dlPrice", i_header = "Region",
#'                            j_header = "LPA_Type", t_header = "Date")
#' summary(est_multi)
#'@export
summary.multi_result <- function(object, ...) {

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
  Y_list <- object$Y_list
  ic <- object$ic
  block_names <- object$block_names
  relative <- matrix(0, 3, R + 1)

  for (i in 1:(R + 1)) {
    if (i < (R + 1)) {
      resid <- rep(1, Ni[i])

      if (r0 == 0) {
        relative[1, i] <- 0
      } else {
        relative[1, i] <- mean(diag(stats::cov(G %*% t(Gamma[[i]])) / (stats::cov(Y_list[[i]]))))
        resid <- resid - diag(stats::cov(G %*% t(Gamma[[i]])) / (stats::cov(Y_list[[i]])))
      }

      if (ri[i] == 0) {
        relative[2, i] <- 0
      } else {
        relative[2, i] <- mean(diag(stats::cov(F[[i]] %*% t(Lambda[[i]])) / (stats::cov(Y_list[[i]]))))
        resid <- resid - diag(stats::cov(F[[i]] %*% t(Lambda[[i]])) / (stats::cov(Y_list[[i]])))
      }

      relative[3, i] <- mean(resid)
    } else {
      relative[, i] <- rowMeans(relative[, 1:R])
    }
  }
  relative<-round(relative, 4)
  sample_size <- c(Ni, N)
  num_local <- c(ri, NA)
  results <- t(rbind(sample_size, num_local, relative))
  rownames(results) <- stringr::str_trunc(c(block_names, "Average/Sum"), 20)
  colnames(results) <- c("Ni", "ri", "RI global", "RI local", "RI error")

  cat(
    "===========================================================",
    "\n"
  )
  cat("\n")
  cat(
    "                  Estimation  results                      ",
    "\n"
  )
  cat("\n")
  cat(
    "===========================================================",
    "\n"
  )
  cat("R = ", R, ", N = ", N, ", T = ", T, " \n")
  cat(r0," global factor(s) detected. r0 = ", r0, " \n")
  cat("local factors are estimated using ", ic, " \n")
  cat("RI = relative importance ratio", "\n")
  print(results)
  cat(
    "===========================================================",
    "\n"
  )
  return(invisible(results))
}



