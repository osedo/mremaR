#'
#' Simulate cell type specific expression analysis
#'
#'This function is based on the simulations described in a benchmarking paper by Meng et al., 2023 (https://doi.org/10.1093/bib/bbac516)
#' @param g number of genes - defaults to 1,000
#' @param n number of samples - defaults to 200, an event spliut between cases and controls
#' @param p proportion of genes that have non-zero LFC in each cell type - defaults to 0.05
#' @param m_lfc mean of the normal distribution from which the absolute lfc values are drawn - dafaults to 4
#' @param sd_lfc sd of the normal distribution from which the absolute lfc values are drawn - dafaults to 0.2
#' @param parameters additional parameters used to simulate expression. See the simulation.parameters document for more details.
#'
#' @returns A list with the following objects:
#'  1) bulk.expression - a g X n matrix of the expression values in mixted samples
#'  2) ols.mixture.estimates - a n X k matrix of cell type proportions estimated using ordinary least squares
#'  3) true.mixture - a n X k matrix of true cell type proportions.
#'  4) ct.reference - a g X k matrix of reference values for cell type specific expression
#'  5) simulated.lfc - a g X k matrix of simulated lfc in each cell type.
#'


CTexpSimulation <- function(g = 1000, n = 200, p = 0.05, m_lfc = 4, sd_lfc = 0.2, parameters = simulation.parameters){

  # number of cell types
  k <- length(parameters$alpha_cases)
  # get lfc values - equal number of cases and controls
  lfc_up_cases <-    stats::rnorm((g*k*p)/2, m_lfc, sd_lfc)
  lfc_up_controls <- stats::rnorm((g*k*p)/2, m_lfc, sd_lfc)
  lfc <- c(lfc_up_cases, -lfc_up_controls)
  # sample mean and dispersion parameters from different cell types
  m_gxk <- MASS::mvrnorm(g, parameters$mean_hat, parameters$mean_covar)
  d_gxk <- MASS::mvrnorm(g, parameters$disp_hat, parameters$disp_covar)

  # cases means
  m_cases_gxk <- m_gxk
  # get random p% of genes from each column
  #select <- unlist(lapply(c(1:k), function(x) g*rep(x-1, g*p))) + unlist(lapply(c(1:k), function(x) sample.int(g, g*p)))
  select <- matrix(data = c(unlist(lapply(c(1:k), function(x) sample.int(g, (g*p)/2))),
                            unlist(lapply(c(1:k), function(x) rep(x, g*p*.5)))), ncol = 2)
  # alter means
  m_cases_gxk[select] <- log((exp(m_cases_gxk[select])) * (2^lfc_up_cases))
  up <- lapply(c(1:k), function(x) c(1:g)[!c(1:g) %in% select[select[,2] == x, 1]])
  select_controls <- matrix(data = c(unlist(lapply(c(1:k), function(x) sample(up[[x]], (g*p)/2))),
                                     unlist(lapply(c(1:k), function(x) rep(x, g*p*.5)))), ncol = 2)
  # alter means
  m_gxk[select_controls] <- log((exp(m_gxk[select_controls])) * (2^lfc_up_controls))
  select <- rbind(select, select_controls)
  #select <- cbind(select, lfc)

  simulated.lfc <- matrix(data = 0, nrow = g, ncol = k)
  simulated.lfc[select] <- lfc

  # simulate the cell type-specific reference panel
  x_controls_gxk <- matrix(stats::rgamma(k*g,
                                  shape = 1/exp(d_gxk),
                                  scale = (exp(d_gxk)*exp(m_gxk))), nrow = g)
  # cases, different mean
  x_cases_gxk <- matrix(stats::rgamma(k*g,
                               shape = 1/exp(d_gxk),
                               scale = (exp(d_gxk)*exp(m_cases_gxk))), nrow = g)

  # sample cell type proportions
  o_cases_k <-    DirichletReg::rdirichlet(n/2, alpha = parameters$alpha_cases)
  o_controls_k <- DirichletReg::rdirichlet(n/2, alpha = parameters$alpha_controls)
  o_k <- rbind(o_controls_k, o_cases_k)
  # get mean expression # weighted sum
  w_controls_gxn <- x_controls_gxk  %*% t(o_controls_k)
  w_cases_gxn <- x_cases_gxk  %*% t(o_cases_k)
  w_gxn <- cbind(w_controls_gxn, w_cases_gxn)
  # simulate observed expression values
  y_gxn <- matrix(stats::rpois(n*g, w_gxn), nrow = g)
  rm(list=setdiff(ls(), c("simulated.lfc", "num", "y_gxn", "o_k", "x_controls_gxk", "x_cases_gxk","select", "g", "k", "n", "p", "m_lfc", "sd_lfc")))

  o_hat_k = apply(y_gxn,2,function(x) stats::lm(x ~ as.matrix(x_controls_gxk))$coefficients[-1])
  o_hat_k = apply(o_hat_k,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
  o_hat_k = apply(o_hat_k,2,function(x) x/sum(x)) #explicit STO constraint
  o_hat_k <- t(o_hat_k)
  decon_rmse <- sqrt(sum(((o_k - o_hat_k)^2)/(n*k)))


  rownames(o_hat_k) <- paste0("sample.", c(1:n))
  colnames(o_hat_k) <- paste0("cellType.", c(1:k))
  rownames(y_gxn) <- paste0("gene.", c(1:g))
  colnames(y_gxn) <- paste0("sample.", c(1:n))
  rownames(o_k) <- paste0("sample.", c(1:n))
  colnames(o_k) <- paste0("cellType.", c(1:k))
  cov <- data.frame("trait" = as.factor(c(rep(0, n/2), rep(1,n/2))))
  row.names(cov) <- paste0("sample.", c(1:n))

  list("bulk.expression" = y_gxn, "ols.mixture.estimate" = o_hat_k, "true.mixture" = o_k,
       "ct.reference" = x_controls_gxk,  "simulated.lfc" = simulated.lfc)


}


#'#' This is data to be included in my package
#'
#' @name data-simulation.parameters
#' @docType data
#' @description
#' A list of default parameters for cell type-specific expression simulation. There are g genes, k cell types, and n samples evenly divided between cases and controls.
#' mean_hat is a vector length k, with mean expression parameters for each cell-types, disp_hat is also a vector of length k with the dispersion parameters for each cell type.
#' mean_covar is a k x k matrix giving the covariance between the mean parameters in each cell type, disp_hat gives the covariance between the dispersion parameters in the different cell types.
#' alpha_cases and alpha_controls are k length vectors with parameters for the dirichlet distribution from which cell type proportion are drawn.
#' These parameters were obtained from the Supplementary materials of Meng et al., 2023 (https://doi.org/10.1093/bib/bbac516)
#'
#' @keywords data
"simulation.parameters"


