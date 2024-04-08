#' Run gene set analysis using Mixture of Random-Effects Meta-Analysis
#'
#' @param postdata A three column dataframe with columns corresponding to gene name, lfc and lfc variance, respectively
#' @param raw.gs A list of named gene sets.
#' @param DF Degrees of freedom between the null and alternative hypotheses.
#' @param threshold The foldchange above which a gene is considered differentially expressed.
#' @param ncores The number of cores to use for parallel computing. Defaults to 1.
#' @param overlap The proportion of the DE normal distributions. Defaults to 0.25
#'
#' @returns A list of two dataframes, results containing the results of the tests and parameters containing the parameters of the distribution of the gene set log fold changes.
#'
#'
#'


mrema <- function(postdata, raw.gs, DF = NULL, threshold = NULL, ncores = 1, overlap = 0.25) {
  # remove genes with na values
  postdata <- postdata[stats::complete.cases(postdata), ]
  # get effect sizes and variance
  effect <- dplyr::pull(postdata, 2)
  variance <- dplyr::pull(postdata, 3)

  # remove gene sets with very low number of genes in dataset
  raw.gs <- raw.gs[unlist(lapply(raw.gs, function(i){
    sum(i %in% dplyr::pull(postdata, 1))
  })) >= 5]
  if(length(raw.gs) < 1) stop("No gene sets with more than 5 genes present in data.")

  ## threshold max for middle component
  comp1_var_max <- seq(0, 5, by = 0.00001)
  comp1_var_max <- comp1_var_max[which(stats::pnorm(log2(threshold), 0, sqrt(comp1_var_max)) < 0.975)[1]]

  # get initial estimates of params from data
  alpha_3 <- sum(stats::pnorm(-log2(threshold), effect, sqrt(variance), lower.tail = TRUE))/nrow(postdata)
  alpha_2 <- sum(stats::pnorm(log2(threshold), effect, sqrt(variance), lower.tail = FALSE))/nrow(postdata)
  alpha_1 <- 1 - (alpha_2 + alpha_3)
  mu2 <- mean(effect[effect >= log2(threshold)])
  mu3 <- mean(effect[effect <= -log2(threshold)])
  mu2 <- ifelse(is.nan(mu2), log2(threshold), mu2)
  mu3 <- ifelse(is.nan(mu3), -log2(threshold), mu3)
  # fit ggm to all genes without regard for set membership
  starting.params <- list("param" = list("mu" = c(0, mu2, mu3), "var" = c(comp1_var_max, 0.5, 0.5), "alpha" = c(alpha_1, alpha_2, alpha_3)))
   print(starting.params)
  all_genes_mixture <- .EM_6FP_fixed(effect, variance, comp1_var_max = comp1_var_max, threshold = threshold, overlap = overlap, starting = starting.params)
  print(all_genes_mixture)
  loglike_all_genes <- all_genes_mixture$loglike


  if (ncores == 1) {
    res <- lapply(seq_along(raw.gs), function(j) .run.mrema(j, raw.gs, postdata, DF, comp1_var_max, threshold, overlap, all_genes_mixture, loglike_all_genes))
  } else {
    clust <- parallel::makeCluster(ncores)
    parallel::clusterEvalQ(clust, {
      library(dplyr)
      library(magrittr)
    })
    parallel::clusterExport(clust, varlist = c(
      "raw.gs", "postdata", "DF", "comp1_var_max", "threshold", "overlap", "all_genes_mixture", "loglike_all_genes", "lower_bound", "tryCatch",
      ".EM_4FP_fixed", ".e_step_set_iter_means", ".m_step_set_iter_means", ".m_step_set_iter_fixed_2DF", ".m_step_set_iter_fixed", ".e_step_set_iter", ".run.mrema", ".EM_6FP_fixed", ".EM_1FP_fixed", ".EM_2FP_fixed", ".e_step_iter", ".m_step_iter_fixed"
    ), envir = environment())
    res <- parallel::parLapply(clust, seq_along(raw.gs), function(j) .run.mrema(j, raw.gs, postdata, DF, comp1_var_max, threshold, overlap, all_genes_mixture, loglike_all_genes))
    parallel::stopCluster(clust)
  }

  analysis.res <- unlist(res, recursive = FALSE)
  params <- unlist(analysis.res[which(names(analysis.res) == "parameters")])

  if (DF == 1) {
    parameters <- data.frame(
      "Gene.Set.mu1" = params[which(names(params) == "parameters.mu1")],
      "Gene.Set.mu2" = params[which(names(params) == "parameters.mu2")],
      "Gene.Set.mu3" = params[which(names(params) == "parameters.mu3")],
      "Gene.Set.var1" = params[which(names(params) == "parameters.var1")],
      "Gene.Set.var2" = params[which(names(params) == "parameters.var2")],
      "Gene.Set.var3" = params[which(names(params) == "parameters.var3")],
      "Gene.Set.alpha1" = params[which(names(params) == "parameters.alpha1")],
      "Gene.Set.alpha2" = params[which(names(params) == "parameters.alpha2")],
      "Gene.Set.alpha3" = params[which(names(params) == "parameters.alpha3")],
      "Backround.mu1" = params[which(names(params) == "parameters.mu1")],
      "Backround.mu2" = params[which(names(params) == "parameters.mu2")],
      "Backround.mu3" = params[which(names(params) == "parameters.mu3")],
      "Backround.var1" = params[which(names(params) == "parameters.var1")],
      "Backround.var2" = params[which(names(params) == "parameters.var2")],
      "Backround.var3" = params[which(names(params) == "parameters.var3")],
      "Backround.alpha1" = params[which(names(params) == "parameters.alpha4")],
      "Backround.alpha2" = params[which(names(params) == "parameters.alpha5")],
      "Backround.alpha3" = params[which(names(params) == "parameters.alpha6")]
    )
  } else if (DF == 2) {
    parameters <- data.frame(
      "Gene.Set.mu1" = params[which(names(params) == "parameters.mu1")],
      "Gene.Set.mu2" = params[which(names(params) == "parameters.mu2")],
      "Gene.Set.mu3" = params[which(names(params) == "parameters.mu3")],
      "Gene.Set.var1" = params[which(names(params) == "parameters.var1")],
      "Gene.Set.var2" = params[which(names(params) == "parameters.var2")],
      "Gene.Set.var3" = params[which(names(params) == "parameters.var3")],
      "Gene.Set.alpha1" = params[which(names(params) == "parameters.alpha1")],
      "Gene.Set.alpha2" = params[which(names(params) == "parameters.alpha2")],
      "Gene.Set.alpha3" = params[which(names(params) == "parameters.alpha3")],
      "Backround.mu1" = params[which(names(params) == "parameters.mu1")],
      "Backround.mu2" = params[which(names(params) == "parameters.mu2")],
      "Backround.mu3" = params[which(names(params) == "parameters.mu3")],
      "Backround.var1" = params[which(names(params) == "parameters.var1")],
      "Backround.var2" = params[which(names(params) == "parameters.var2")],
      "Backround.var3" = params[which(names(params) == "parameters.var3")],
      "Backround.alpha1" = params[which(names(params) == "parameters.alpha4")],
      "Backround.alpha2" = params[which(names(params) == "parameters.alpha5")],
      "Backround.alpha3" = params[which(names(params) == "parameters.alpha6")]
    )
  } else if (DF == 6) {
    parameters <- data.frame(
      "Gene.Set.mu1" = params[which(names(params) == "parameters.Gene.Set.mu1")],
      "Gene.Set.mu2" = params[which(names(params) == "parameters.Gene.Set.mu2")],
      "Gene.Set.mu3" = params[which(names(params) == "parameters.Gene.Set.mu3")],
      "Gene.Set.var1" = params[which(names(params) == "parameters.Gene.Set.var1")],
      "Gene.Set.var2" = params[which(names(params) == "parameters.Gene.Set.var2")],
      "Gene.Set.var3" = params[which(names(params) == "parameters.Gene.Set.var3")],
      "Gene.Set.alpha1" = params[which(names(params) == "parameters.Gene.Set.alpha1")],
      "Gene.Set.alpha2" = params[which(names(params) == "parameters.Gene.Set.alpha2")],
      "Gene.Set.alpha3" = params[which(names(params) == "parameters.Gene.Set.alpha3")],
      "Backround.mu1" = params[which(names(params) == "parameters.Backround.mu1")],
      "Backround.mu2" = params[which(names(params) == "parameters.Backround.mu2")],
      "Backround.mu3" = params[which(names(params) == "parameters.Backround.mu3")],
      "Backround.var1" = params[which(names(params) == "parameters.Backround.var1")],
      "Backround.var2" = params[which(names(params) == "parameters.Backround.var2")],
      "Backround.var3" = params[which(names(params) == "parameters.Backround.var3")],
      "Backround.alpha1" = params[which(names(params) == "parameters.Backround.alpha1")],
      "Backround.alpha2" = params[which(names(params) == "parameters.Backround.alpha2")],
      "Backround.alpha3" = params[which(names(params) == "parameters.Backround.alpha3")]
    )
  } else if (DF == 4) {
    parameters <- data.frame(
      "Gene.Set.mu1" = params[which(names(params) == "parameters.mu1")],
      "Gene.Set.mu2" = params[which(names(params) == "parameters.mu2")],
      "Gene.Set.mu3" = params[which(names(params) == "parameters.mu3")],
      "Gene.Set.var1" = params[which(names(params) == "parameters.var1")],
      "Gene.Set.var2" = params[which(names(params) == "parameters.var2")],
      "Gene.Set.var3" = params[which(names(params) == "parameters.var3")],
      "Gene.Set.alpha1" = params[which(names(params) == "parameters.alpha1")],
      "Gene.Set.alpha2" = params[which(names(params) == "parameters.alpha2")],
      "Gene.Set.alpha3" = params[which(names(params) == "parameters.alpha3")],
      "Backround.mu1" = params[which(names(params) == "parameters.mu4")],
      "Backround.mu2" = params[which(names(params) == "parameters.mu5")],
      "Backround.mu3" = params[which(names(params) == "parameters.mu6")],
      "Backround.var1" = params[which(names(params) == "parameters.var1")],
      "Backround.var2" = params[which(names(params) == "parameters.var2")],
      "Backround.var3" = params[which(names(params) == "parameters.var3")],
      "Backround.alpha1" = params[which(names(params) == "parameters.alpha4")],
      "Backround.alpha2" = params[which(names(params) == "parameters.alpha5")],
      "Backround.alpha3" = params[which(names(params) == "parameters.alpha6")]
    )
  }
  res <- do.call("rbind", analysis.res[which(names(analysis.res) == "test")])
  res$ADJ.PVAL <- stats::p.adjust(res$PVAL, method = "BH", n = nrow(res))
  rownames(res) <- NULL
  detected_gene_sets <- res
  rownames(parameters) <- res$GENE.SET
  list("results" = detected_gene_sets, "parameters" = parameters)

  # if(is.null(params) == TRUE) return(detected_gene_sets) else return(parameters_h1)
}


.run.mrema <- function(j, raw.gs, postdata, DF, comp1_var_max, threshold, overlap, all_genes_mixture, loglike_all_genes) {
  tryCatch(
    {
      if (DF == 1) {
        ### run 1DF test
        set_specific_post <- postdata
        set_specific_post$set <- ifelse(dplyr::pull(postdata, 1) %in% raw.gs[[j]], 1, 0)
        effect_set <- dplyr::pull(set_specific_post, 2)
        variance_set <- dplyr::pull(set_specific_post, 3)
        set <- set_specific_post$set
        set_mixture <- .EM_1FP_fixed(effect_set, variance_set, set, comp1_var_max, threshold = threshold, overlap = overlap, starting = all_genes_mixture)
        loglike_set_genes <- set_mixture$loglike
        set_parameters <- set_mixture$param
        ll_trace <- set_mixture$ll.vector
        # compare the two models
        teststat <- 2*(-loglike_all_genes - (-loglike_set_genes))
        pval <- stats::pchisq(teststat,df=1,lower.tail=FALSE)

        BIC_all <- 6 * log(nrow(postdata)) - 2 * loglike_all_genes
        BIC_set <- 7 * log(nrow(postdata)) - 2 * (loglike_set_genes)
        parameters_h1 <- set_parameters
        nonDE_criterion <- set_parameters$alpha[1] < set_parameters$alpha[4]
        weight_diff <- (set_parameters$alpha[4] - set_parameters$alpha[1])
        v <- pval
        tol <- length(raw.gs[[j]])
        BIC <- BIC_set < BIC_all
        # progress(j, max.value = length(raw.gs))
      } else if (DF == 6) {
        ### run 6DF approach
        set_specific_post <- postdata
        set_specific_post$set <- ifelse(dplyr::pull(postdata, 1) %in% raw.gs[[j]], 1, 0)
        set_specific_post_in <- set_specific_post %>% dplyr::filter(set == 1)
        effect_inset <- dplyr::pull(set_specific_post_in, 2)
        variance_inset <- dplyr::pull(set_specific_post_in, 3)

        # fit the gmm to the genes in the gene set
        # get initial estimates of params from data
        alpha_3 <- sum(stats::pnorm(-log2(threshold), effect_inset, sqrt(variance_inset), lower.tail = TRUE))/length(effect_inset)
        alpha_2 <- sum(stats::pnorm(log2(threshold), effect_inset, sqrt(variance_inset), lower.tail = FALSE))/length(effect_inset)
        alpha_1 <- 1 - (alpha_2 + alpha_3)
        mu2 <- mean(effect_inset[effect_inset >= log2(threshold)])
        mu3 <- mean(effect_inset[effect_inset <= -log2(threshold)])
        mu2 <- ifelse(is.nan(mu2), log2(threshold), mu2)
        mu3 <- ifelse(is.nan(mu3), -log2(threshold), mu3)
        # fit ggm to all genes without regard for set membership
        starting.params <- list("param" = list("mu" = c(0, mu2, mu3), "var" = c(comp1_var_max, 0.5, 0.5), "alpha" = c(alpha_1, alpha_2, alpha_3)))
        inset_mixture <- .EM_6FP_fixed(effect_inset, variance_inset, comp1_var_max = comp1_var_max, threshold = threshold, overlap = overlap, starting = starting.params)
        loglike_Inset_genes <- inset_mixture$loglike
        inset_parameters <- inset_mixture$param

        # get all genes in the outset
        set_specific_post_out <- set_specific_post %>% dplyr::filter(set == 0)
        effect_outset <- dplyr::pull(set_specific_post_out, 2)
        variance_outset <- dplyr::pull(set_specific_post_out, 3)

        # fit the gmm to genes outside the gene set
        starting.params <- list("param" = list("mu" = c(0, log2(threshold) + 1, -log2(threshold) - 1), "var" = c(0.5, 0.5, 0.5), "alpha" = c(0.8, 0.1, 0.1)))
        outset_mixture <- .EM_6FP_fixed(effect_outset, variance_outset, comp1_var_max = comp1_var_max, threshold = threshold, overlap = overlap, starting = starting.params)
        loglike_Outset_genes <- outset_mixture$loglike
        outset_parameters <- outset_mixture$param
        parameters_h1 <- list("Gene.Set" = inset_parameters, "Backround" = outset_parameters)

        # compare the two models
        teststat <- 2*(-loglike_all_genes - (-loglike_Inset_genes-loglike_Outset_genes))
        pval <- stats::pchisq(teststat,df=6,lower.tail=FALSE)
        BIC_all <- 6 * log(nrow(postdata)) - 2 * loglike_all_genes
        BIC_set <- 12 * log(nrow(postdata)) - 2 * (loglike_Inset_genes + loglike_Outset_genes)
        v <- pval
       # print(pval)
        nonDE_criterion <- inset_parameters$alpha[1] < outset_parameters$alpha[1]
        weight_diff <- (outset_parameters$alpha[1] - inset_parameters$alpha[1])
        tol <- length(raw.gs[[j]])
        BIC <- BIC_set < BIC_all

        #  progress(j, max.value = length(raw.gs))
      } else if (DF == 2) {
        ### run 2DF approach
        set_specific_post <- postdata
        set_specific_post$set <- ifelse(dplyr::pull(postdata, 1) %in% raw.gs[[j]], 1, 0)
        effect_set <- dplyr::pull(set_specific_post, 2)
        variance_set <- dplyr::pull(set_specific_post, 3)
        set <- set_specific_post$set
        set_mixture <- .EM_2FP_fixed(effect_set, variance_set, set, comp1_var_max, threshold = threshold, overlap = overlap, starting = all_genes_mixture)
        loglike_set_genes <- set_mixture$loglike
        set_parameters <- set_mixture$param
        ll_trace <- set_mixture$ll.vector
        # compare the two models
        teststat <- 2*(-loglike_all_genes - (-loglike_set_genes))
        pval <- stats::pchisq(teststat,df=2,lower.tail=FALSE)
        BIC_all <- 6 * log(nrow(postdata)) - 2 * loglike_all_genes
        BIC_set <- 8 * log(nrow(postdata)) - 2 * (loglike_set_genes)
        parameters_h1 <- set_parameters
        nonDE_criterion <- set_parameters$alpha[1] < set_parameters$alpha[4]
        weight_diff <- (set_parameters$alpha[4] - set_parameters$alpha[1])
        v <- pval
        tol <- length(raw.gs[[j]])
        BIC <- BIC_set < BIC_all
        # progress(j, max.value = length(raw.gs))
      } else if (DF == 4) {
        # print("Running 4 DF")
        ### run 2DF approach
        set_specific_post <- postdata
        set_specific_post$set <- ifelse(dplyr::pull(postdata, 1) %in% raw.gs[[j]], 1, 0)
        effect_set <- dplyr::pull(set_specific_post, 2)
        variance_set <- dplyr::pull(set_specific_post, 3)
        set <- set_specific_post$set
        set_mixture <- .EM_4FP_fixed(effect_set, variance_set, set, comp1_var_max, threshold = threshold, overlap = overlap, starting = all_genes_mixture)
        # print(set_mixture)
        loglike_set_genes <- set_mixture$loglike
        set_parameters <- set_mixture$param
        ll_trace <- set_mixture$ll.vector
        # compare the two models
        teststat <- 2*(-loglike_all_genes - (-loglike_set_genes))
        pval <- stats::pchisq(teststat,df=4,lower.tail=FALSE)
        BIC_all <- 6 * log(nrow(postdata)) - 2 * loglike_all_genes
        BIC_set <- 10 * log(nrow(postdata)) - 2 * (loglike_set_genes)
        parameters_h1 <- set_parameters

        nonDE_criterion <- set_parameters$alpha[1] < set_parameters$alpha[4]
        weight_diff <- (set_parameters$alpha[4] - set_parameters$alpha[1])
        v <- pval
        tol <- length(raw.gs[[j]])
        BIC <- BIC_set < BIC_all
      }
    },
    error = function(e) {
      cat("Error in ", names(raw.gs[j]), "set\n")
      warning(e)
      nonDE_criterion <- NA
      v <- NA
      tol <- NA
      BIC <- NA
      # progress(j, max.value = length(raw.gs))
    }
  )
  res <- data.frame("GENE.SET" = names(raw.gs)[j], "Prop.DE.Increased" = as.numeric(nonDE_criterion), "Estimated.Difference" = weight_diff, "NR.GENES" = tol, "PVAL" = v)
  list("test" = res, "parameters" = parameters_h1)
}


# fitting the ggm with all parameters (apart from non-DE component mean and variance) free
.EM_6FP_fixed <- function(effect, variance, effect_summary_df, comp1_var_max, threshold, overlap = overlap, starting) {
  n <- 1000
  m <- 1e-6
  iter <- 10
  for (i in 1:n) {
    if (i == 1) {
      # Initialization
      e.step <- .e_step_iter(effect, variance, starting$param$mu, starting$param$var, c(starting$param$alpha))
      m.step <- .m_step_iter_fixed(effect, variance, starting$param$var, iter, e.step[["posterior_df"]], comp1_var_max, threshold, overlap = overlap)
      print(m.step)
      cur.loglik <- e.step[["loglik"]]
      loglik.vector <- e.step[["loglik"]]
    } else {
      # Repeat E and M steps till convergence
      e.step <- .e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], m.step[["alpha"]])
      m.step <- .m_step_iter_fixed(effect, variance, m.step[["var"]], iter, e.step[["posterior_df"]], comp1_var_max, threshold, overlap = overlap)
      print(m.step)
      loglik.vector <- c(loglik.vector, e.step[["loglik"]])
      loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
      if (loglik.diff < m) {
        break
      } else {
        cur.loglik <- e.step[["loglik"]]
      }
    }
  }
  loglike_all_genes <- utils::tail(loglik.vector, n = 1)
  parameters <- list("loglike" = loglike_all_genes, "param" = m.step)
  return(parameters)
}

## # fitting the ggm 1DF approach
.EM_1FP_fixed <- function(effect_set, variance_set, set, comp1_var_max, threshold, overlap = overlap, starting) {
  n <- 1000
  m <- 1e-6
  iter <- 10
  for (i in 1:n) {
    if (i == 1) {
      # Initialization
      e.step.set <- .e_step_set_iter(effect_set, variance_set, set, starting$param$mu, starting$param$var, c(starting$param$alpha, starting$param$alpha))

      m.step.set <- .m_step_set_iter_fixed(effect_set, variance_set, set, starting$param$var, iter, e.step.set[["posterior_df"]], comp1_var_max, threshold = threshold, overlap = overlap)
      cur.loglik.set <- e.step.set[["loglik"]]
      loglik.vector.set <- e.step.set[["loglik"]]
    } else {
      # Repeat E and M steps till convergence
      e.step.set <- .e_step_set_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], m.step.set[["alpha"]])
      m.step.set <- .m_step_set_iter_fixed(effect_set, variance_set, set, m.step.set[["var"]], iter, e.step.set[["posterior_df"]], comp1_var_max, threshold = threshold, overlap = overlap)

      loglik.vector.set <- c(loglik.vector.set, e.step.set[["loglik"]])

      loglik.diff.set <- abs((cur.loglik.set - e.step.set[["loglik"]]))

      if (loglik.diff.set < m) {
        break
      } else {
        cur.loglik.set <- e.step.set[["loglik"]]
      }
    }
  }

  loglike_set_genes <- utils::tail(loglik.vector.set, n = 1)
  parameters <- list("loglike" = loglike_set_genes, "param" = m.step.set, "ll.vector" = loglik.vector.set)
  return(parameters)
}

## # fitting the ggm 1DF approach
.EM_4FP_fixed <- function(effect_set, variance_set, set, comp1_var_max, threshold, overlap = overlap, starting) {
  n <- 1000
  m <- 1e-6
  iter <- 10
  for (i in 1:n) {
    if (i == 1) {
      # Initialization
      e.step.set <- .e_step_set_iter_means(effect_set, variance_set, set, c(starting$param$mu, starting$param$mu), starting$param$var, c(starting$param$alpha, starting$param$alpha))

      m.step.set <- .m_step_set_iter_means(effect_set, variance_set, set, starting$param$var, iter, e.step.set[["posterior_df"]], comp1_var_max, threshold = threshold, overlap = overlap)
      # print(m.step.set)
      cur.loglik.set <- e.step.set[["loglik"]]
      loglik.vector.set <- e.step.set[["loglik"]]
    } else {
      # Repeat E and M steps till convergence
      e.step.set <- .e_step_set_iter_means(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], m.step.set[["alpha"]])
      m.step.set <- .m_step_set_iter_means(effect_set, variance_set, set, m.step.set[["var"]], iter, e.step.set[["posterior_df"]], comp1_var_max, threshold = threshold, overlap = overlap)

      loglik.vector.set <- c(loglik.vector.set, e.step.set[["loglik"]])

      loglik.diff.set <- abs((cur.loglik.set - e.step.set[["loglik"]]))

      if (loglik.diff.set < m) {
        break
      } else {
        cur.loglik.set <- e.step.set[["loglik"]]
      }
    }
  }

  loglike_set_genes <- utils::tail(loglik.vector.set, n = 1)
  parameters <- list("loglike" = loglike_set_genes, "param" = m.step.set, "ll.vector" = loglik.vector.set)
  return(parameters)
}




# the e_step and m_step functions in the EM algorithm
.e_step_iter <- function(x, lfc_var, mu_vector, component_var, alpha_vector) {
  # both lfc_se and component_var contribute to total variance

  comp1_prod <- stats::dnorm(x, mu_vector[1], sqrt(component_var[1] + lfc_var)) * alpha_vector[1]
  comp2_prod <- stats::dnorm(x, mu_vector[2], sqrt(component_var[2] + lfc_var)) * alpha_vector[2]
  comp3_prod <- stats::dnorm(x, mu_vector[3], sqrt(component_var[3] + lfc_var)) * alpha_vector[3]

  sum_of_comps <- comp1_prod + comp2_prod + comp3_prod
  sum_of_comps[which(sum_of_comps == 0)] <- 1e-200
  comp1_post <- comp1_prod / sum_of_comps
  comp2_post <- comp2_prod / sum_of_comps
  comp3_post <- comp3_prod / sum_of_comps

  sum_of_comps_ln <- log(sum_of_comps, base = exp(1))
  sum_of_comps_ln_sum <- sum(sum_of_comps_ln)

  list(
    "loglik" = sum_of_comps_ln_sum,
    "posterior_df" = cbind(comp1_post, comp2_post, comp3_post),
    "prod_df" = cbind(comp1_prod, comp2_prod, comp3_prod)
  )
}


# the e_step and m_step functions in the EM algorithm for w_o_mrema (weights only)
.e_step_set_iter <- function(x, lfc_var, set, mu_vector, component_var, alpha_vector) {
  # both sd_vector and sigma are variance


  comp1_prod <- stats::dnorm(x, mu_vector[1], sqrt(component_var[1] + lfc_var)) * alpha_vector[1] * set
  comp2_prod <- stats::dnorm(x, mu_vector[2], sqrt(component_var[2] + lfc_var)) * alpha_vector[2] * set
  comp3_prod <- stats::dnorm(x, mu_vector[3], sqrt(component_var[3] + lfc_var)) * alpha_vector[3] * set
  # posteiror dist for genes outside of the set

  comp4_prod <- stats::dnorm(x, mu_vector[1], sqrt(component_var[1] + lfc_var)) * alpha_vector[4] * (1 - set)
  comp5_prod <- stats::dnorm(x, mu_vector[2], sqrt(component_var[2] + lfc_var)) * alpha_vector[5] * (1 - set)
  comp6_prod <- stats::dnorm(x, mu_vector[3], sqrt(component_var[3] + lfc_var)) * alpha_vector[6] * (1 - set)

  sum_of_comps1 <- comp1_prod + comp2_prod + comp3_prod
  sum_of_comps2 <- comp4_prod + comp5_prod + comp6_prod
  sum_of_comps <- sum_of_comps1 + sum_of_comps2
  sum_of_comps[which(sum_of_comps == 0)] <- 1e-200

  comp1_post <- comp1_prod / sum_of_comps1
  comp1_post[is.na(comp1_post)] <-  1e-200

  comp2_post <- comp2_prod / sum_of_comps1
  comp2_post[is.na(comp2_post)] <-  1e-200

  comp3_post <- comp3_prod / sum_of_comps1
  comp3_post[is.na(comp3_post)] <-  1e-200

  comp4_post <- comp4_prod / sum_of_comps2
  comp4_post[is.na(comp4_post)] <-  1e-200

  comp5_post <- comp5_prod / sum_of_comps2
  comp5_post[is.na(comp5_post)] <-  1e-200

  comp6_post <- comp6_prod / sum_of_comps2
  comp6_post[is.na(comp6_post)] <-  1e-200

  sum_of_comps_ln <- log(sum_of_comps, base = exp(1))
  sum_of_comps_ln_sum <- sum(sum_of_comps_ln)

  list(
    "loglik" = sum_of_comps_ln_sum,
    "posterior_df" = cbind(comp1_post, comp2_post, comp3_post, comp4_post, comp5_post, comp6_post)
  )
}




# the e_step and m_step functions in the EM algorithm for w_o_mrema (weights only)
.e_step_set_iter_means <- function(x, lfc_var, set, mu_vector, component_var, alpha_vector) {
  # both sd_vector and sigma are variance


  comp1_prod <- stats::dnorm(x, mu_vector[1], sqrt(component_var[1] + lfc_var)) * alpha_vector[1] * set
  comp2_prod <- stats::dnorm(x, mu_vector[2], sqrt(component_var[2] + lfc_var)) * alpha_vector[2] * set
  comp3_prod <- stats::dnorm(x, mu_vector[3], sqrt(component_var[3] + lfc_var)) * alpha_vector[3] * set
  # posteiror dist for genes outside of the set

  comp4_prod <- stats::dnorm(x, mu_vector[4], sqrt(component_var[1] + lfc_var)) * alpha_vector[4] * (1 - set)
  comp5_prod <- stats::dnorm(x, mu_vector[5], sqrt(component_var[2] + lfc_var)) * alpha_vector[5] * (1 - set)
  comp6_prod <- stats::dnorm(x, mu_vector[6], sqrt(component_var[3] + lfc_var)) * alpha_vector[6] * (1 - set)

  sum_of_comps1 <- comp1_prod + comp2_prod + comp3_prod
  sum_of_comps2 <- comp4_prod + comp5_prod + comp6_prod
  sum_of_comps <- sum_of_comps1 + sum_of_comps2
  sum_of_comps[which(sum_of_comps == 0)] <- 1e-200

  comp1_post <- comp1_prod / sum_of_comps1
  comp1_post[is.na(comp1_post)] <-  1e-200

  comp2_post <- comp2_prod / sum_of_comps1
  comp2_post[is.na(comp2_post)] <-  1e-200

  comp3_post <- comp3_prod / sum_of_comps1
  comp3_post[is.na(comp3_post)] <-  1e-200

  comp4_post <- comp4_prod / sum_of_comps2
  comp4_post[is.na(comp4_post)] <-  1e-200

  comp5_post <- comp5_prod / sum_of_comps2
  comp5_post[is.na(comp5_post)] <-  1e-200

  comp6_post <- comp6_prod / sum_of_comps2
  comp6_post[is.na(comp6_post)] <-  1e-200

  sum_of_comps_ln <- log(sum_of_comps, base = exp(1))
  sum_of_comps_ln_sum <- sum(sum_of_comps_ln)

  list(
    "loglik" = sum_of_comps_ln_sum,
    "posterior_df" = cbind(comp1_post, comp2_post, comp3_post, comp4_post, comp5_post, comp6_post)
  )
}


.m_step_set_iter_means <- function(x, lfc_var, set, component_var, t, posterior_df, comp1_var_max, threshold, overlap = overlap) {
  comp1_n <- sum(posterior_df[, 1])
  comp2_n <- sum(posterior_df[, 2])
  comp3_n <- sum(posterior_df[, 3])
  comp4_n <- sum(posterior_df[, 4])
  comp5_n <- sum(posterior_df[, 5])
  comp6_n <- sum(posterior_df[, 6])

  #  comp1_4_n<- comp1_n + comp4_n
  #  comp2_5_n<- comp2_n + comp5_n
  #  comp3_6_n<- comp3_n + comp6_n
  ###########################
  for (i in 1:t) {
    if (i == 1) {
      comp2_var <- component_var[2]
      w2_i <- 1 / (comp2_var + lfc_var)

      comp2_mean_min <- stats::qnorm(1 - overlap, log2(threshold), sqrt(comp2_var))
      comp2_mu <- max(comp2_mean_min, sum(posterior_df[, 2] * w2_i * x) / sum(posterior_df[, 2] * w2_i))
      comp5_mu <- max(comp2_mean_min, sum(posterior_df[, 5] * w2_i * x) / sum(posterior_df[, 5] * w2_i))

      comp3_var <- component_var[3]
      w3_i <- 1 / (comp3_var + lfc_var)
      comp3_mean_max <- stats::qnorm(overlap, -log2(threshold), sqrt(comp3_var))
      comp3_mu <- min(comp3_mean_max, sum(posterior_df[, 3] * w3_i * x) / sum(posterior_df[, 3] * w3_i))
      comp6_mu <- min(comp3_mean_max, sum(posterior_df[, 6] * w3_i * x) / sum(posterior_df[, 6] * w3_i))

      comp1_var <- comp1_var_max

      comp1_mu <- 0
      comp4_mu <- 0
    } else {
      comp2_var <- max(0.005, (sum(posterior_df[, 2] * ((w2_i)^2) * (((x - comp2_mu)^2) - lfc_var)) + sum(posterior_df[, 5] * ((w2_i)^2) * (((x - comp5_mu)^2) - lfc_var))) * (1 / sum((posterior_df[, 2] + posterior_df[, 5]) * w2_i^2)))
      w2_i <- 1 / (comp2_var + lfc_var)
      comp2_mean_min <- stats::qnorm((1 - overlap), log2(threshold), sqrt(comp2_var))
      comp2_mu <- max(comp2_mean_min, sum(posterior_df[, 2] * w2_i * x) / sum(posterior_df[, 2] * w2_i))
      comp5_mu <- max(comp2_mean_min, sum(posterior_df[, 5] * w2_i * x) / sum(posterior_df[, 5] * w2_i))

      comp3_var <- max(0.005, (sum(posterior_df[, 3] * ((w3_i)^2) * (((x - comp3_mu)^2) - lfc_var)) + sum(posterior_df[, 6] * ((w3_i)^2) * (((x - comp6_mu)^2) - lfc_var))) * (1 / sum((posterior_df[, 3] + posterior_df[, 6]) * w3_i^2)))
      w3_i <- 1 / (comp3_var + lfc_var)
      comp3_mean_max <- stats::qnorm(overlap, -log2(threshold), sqrt(comp3_var))
      comp3_mu <- min(comp3_mean_max, sum(posterior_df[, 3] * w3_i * x) / sum(posterior_df[, 3] * w3_i))
      comp6_mu <- min(comp3_mean_max, sum(posterior_df[, 6] * w3_i * x) / sum(posterior_df[, 6] * w3_i))

      comp1_var <- comp1_var_max
      comp1_mu <- 0
      comp4_mu <- 0
    }
  }


  comp1_alpha <- max(comp1_n / sum(set), lower_bound)
  comp2_alpha <- max(comp2_n / sum(set), lower_bound)
  comp3_alpha <- max(comp3_n / sum(set), lower_bound)
  comp4_alpha <- max(comp4_n / (length(set) - sum(set)), lower_bound)
  comp5_alpha <- max(comp5_n / (length(set) - sum(set)), lower_bound)
  comp6_alpha <- max(comp6_n / (length(set) - sum(set)), lower_bound)

  list(
    "mu" = c(comp1_mu, comp2_mu, comp3_mu, comp4_mu, comp5_mu, comp6_mu),
    "var" = c(comp1_var, comp2_var, comp3_var),
    "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha, comp6_alpha)
  )
}


.m_step_iter_fixed <- function(x, lfc_var, component_var, t, posterior_df, comp1_var_max, threshold, overlap = overlap) {
  comp1_n <- sum(posterior_df[, 1])
  comp2_n <- sum(posterior_df[, 2])
  comp3_n <- sum(posterior_df[, 3])


  ###########################
  for (i in 1:t) {
    if (i == 1) {
      # Initialization
      comp2_var <- component_var[2]
      ## weights
      w2_i <- 1 / (comp2_var + lfc_var)

      # gets the mininum mean allowed when variance is comp2_var to keep 75% of the distribution above threshold
      comp2_mean_min <- stats::qnorm((1 - overlap), log2(threshold), sqrt(comp2_var))
      ## use either minimum value or mean estimate if bigger
      comp2_mu <- max(comp2_mean_min, sum(posterior_df[, 2] * w2_i * x) / sum(posterior_df[, 2] * w2_i))

      comp3_var <- component_var[3]
      ## weights
      w3_i <- 1 / (comp3_var + lfc_var)
      # gets the maximum mean allowed when variance is comp3_var to keep 75% of the distribution below threshold
      comp3_mean_max <- stats::qnorm(overlap, -log2(threshold), sqrt(comp3_var))
      ## use either maximum value or mean estimate if smaller
      comp3_mu <- min(comp3_mean_max, sum(posterior_df[, 3] * w3_i * x) / sum(posterior_df[, 3] * w3_i))

      comp1_var <- comp1_var_max

      comp1_mu <- 0
    } else {
      comp2_var <- max(0.005, sum(posterior_df[, 2] * ((w2_i)^2) * (((x - comp2_mu)^2) - lfc_var)) / (sum(posterior_df[, 2] * w2_i^2)))
      w2_i <- 1 / (comp2_var + lfc_var)
      comp2_mean_min <- stats::qnorm((1 - overlap), log2(threshold), sqrt(comp2_var))
      comp2_mu <- max(comp2_mean_min, sum(posterior_df[, 2] * w2_i * x) / sum(posterior_df[, 2] * w2_i))

      comp3_var <- max(0.005, sum(posterior_df[, 3] * ((w3_i)^2) * (((x - comp3_mu)^2) - lfc_var)) / (sum(posterior_df[, 3] * w3_i^2)))
      w3_i <- 1 / (comp3_var + lfc_var)
      comp3_mean_max <- stats::qnorm(overlap, -log2(threshold), sqrt(comp3_var))
      comp3_mu <- min(comp3_mean_max, sum(posterior_df[, 3] * w3_i * x) / sum(posterior_df[, 3] * w3_i))

      comp1_var <- comp1_var_max

      comp1_mu <- 0
    }
  }



  comp1_alpha <- max(comp1_n / length(x), lower_bound)
  comp2_alpha <- max(comp2_n / length(x), lower_bound)
  comp3_alpha <- max(comp3_n / length(x), lower_bound)

  list(
    "mu" = c(comp1_mu, comp2_mu, comp3_mu),
    "var" = c(comp1_var, comp2_var, comp3_var),
    "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha)
  )
}


.m_step_set_iter_fixed <- function(x, lfc_var, set, component_var, t, posterior_df, comp1_var_max, threshold, overlap = overlap) {
  comp1_n <- sum(posterior_df[, 1])
  comp2_n <- sum(posterior_df[, 2])
  comp3_n <- sum(posterior_df[, 3])
  comp4_n <- sum(posterior_df[, 4])
  comp5_n <- sum(posterior_df[, 5])
  comp6_n <- sum(posterior_df[, 6])

  comp1_alpha <- max(comp1_n / sum(set), lower_bound)
  comp4_alpha <- max(comp4_n / (length(set) - sum(set)), lower_bound)

  c <- (comp2_n + comp5_n) / (comp2_n + comp5_n + comp3_n + comp6_n)

  comp2_alpha <- (1 - comp1_alpha) * (c)
  comp3_alpha <- (1 - comp1_alpha) * (1 - c)

  comp5_alpha <- (1 - comp4_alpha) * (c)
  comp6_alpha <- (1 - comp4_alpha) * (1 - c)



  # the following for loop iterates between mean and variance, the number of iteration is t,
  comp1_pd <- posterior_df[, 1] + posterior_df[, 4]
  comp2_pd <- posterior_df[, 2] + posterior_df[, 5]
  comp3_pd <- posterior_df[, 3] + posterior_df[, 6]
  ###########################
  for (i in 1:t) {
    if (i == 1) {
      comp2_var <- component_var[2]
      w2_i <- 1 / (comp2_var + lfc_var)
      comp2_mean_min <- stats::qnorm(1 - overlap, log2(threshold), sqrt(comp2_var))
      comp2_mu <- max(comp2_mean_min, sum(comp2_pd * w2_i * x) / sum(comp2_pd * w2_i))

      comp3_var <- component_var[3]
      w3_i <- 1 / (comp3_var + lfc_var)
      comp3_mean_max <- stats::qnorm(overlap, -log2(threshold), sqrt(comp3_var))
      comp3_mu <- min(comp3_mean_max, sum(comp3_pd * w3_i * x) / sum(comp3_pd * w3_i))

      comp1_var <- comp1_var_max

      comp1_mu <- 0
    } else {
      comp2_var <- max(0.005, sum(comp2_pd * ((w2_i)^2) * (((x - comp2_mu)^2) - lfc_var)) / (sum(comp2_pd * w2_i^2)))
      w2_i <- 1 / (comp2_var + lfc_var)
      comp2_mean_min <- stats::qnorm((1 - overlap), log2(threshold), sqrt(comp2_var))
      comp2_mu <- max(comp2_mean_min, sum(comp2_pd * w2_i * x) / sum(comp2_pd * w2_i))

      comp3_var <- max(0.005, sum(comp3_pd * ((w3_i)^2) * (((x - comp3_mu)^2) - lfc_var)) / (sum(comp3_pd * w3_i^2)))
      w3_i <- 1 / (comp3_var + lfc_var)
      comp3_mean_max <- stats::qnorm(overlap, -log2(threshold), sqrt(comp3_var))
      comp3_mu <- min(comp3_mean_max, sum(comp3_pd * w3_i * x) / sum(comp3_pd * w3_i))

      comp1_var <- comp1_var_max
      comp1_mu <- 0
    }
  }


  comp1_alpha <- max(comp1_alpha, lower_bound)
  comp2_alpha <- max(comp2_alpha, lower_bound)
  comp3_alpha <- max(comp3_alpha, lower_bound)
  comp4_alpha <- max(comp4_alpha, lower_bound)
  comp5_alpha <- max(comp5_alpha, lower_bound)
  comp6_alpha <- max(comp6_alpha, lower_bound)

  list(
    "mu" = c(comp1_mu, comp2_mu, comp3_mu),
    "var" = c(comp1_var, comp2_var, comp3_var),
    "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha, comp6_alpha)
  )
}



## # fitting the ggm 2DF approach
.EM_2FP_fixed <- function(effect_set, variance_set, set, comp1_var_max, threshold, overlap = overlap, starting) {
  n <- 1000
  m <- 1e-6
  iter <- 10
  for (i in 1:n) {
    if (i == 1) {
      # Initialization
      e.step.set <- .e_step_set_iter(effect_set, variance_set, set, starting$param$mu, starting$param$var, c(starting$param$alpha, starting$param$alpha))

      m.step.set <- .m_step_set_iter_fixed_2DF(effect_set, variance_set, set, starting$param$var, iter, e.step.set[["posterior_df"]], comp1_var_max, threshold = threshold, overlap = overlap)
      cur.loglik.set <- e.step.set[["loglik"]]
      loglik.vector.set <- e.step.set[["loglik"]]
    } else {
      # Repeat E and M steps till convergence
      e.step.set <- .e_step_set_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], m.step.set[["alpha"]])
      m.step.set <- .m_step_set_iter_fixed_2DF(effect_set, variance_set, set, m.step.set[["var"]], iter, e.step.set[["posterior_df"]], comp1_var_max, threshold = threshold, overlap = overlap)

      loglik.vector.set <- c(loglik.vector.set, e.step.set[["loglik"]])

      loglik.diff.set <- abs((cur.loglik.set - e.step.set[["loglik"]]))

      if (loglik.diff.set < m) {
        break
      } else {
        cur.loglik.set <- e.step.set[["loglik"]]
      }
    }
  }

  loglike_set_genes <- utils::tail(loglik.vector.set, n = 1)
  parameters <- list("loglike" = loglike_set_genes, "param" = m.step.set, "ll.vector" = loglik.vector.set)
  return(parameters)
}


.m_step_set_iter_fixed_2DF <- function(x, lfc_var, set, component_var, t, posterior_df, comp1_var_max, threshold, overlap = overlap) {
  comp1_n <- sum(posterior_df[, 1])
  comp2_n <- sum(posterior_df[, 2])
  comp3_n <- sum(posterior_df[, 3])
  comp4_n <- sum(posterior_df[, 4])
  comp5_n <- sum(posterior_df[, 5])
  comp6_n <- sum(posterior_df[, 6])





  # the following for loop iterates between mean and variance, the number of iteration is t,
  comp1_pd <- posterior_df[, 1] + posterior_df[, 4]
  comp2_pd <- posterior_df[, 2] + posterior_df[, 5]
  comp3_pd <- posterior_df[, 3] + posterior_df[, 6]
  ###########################
  for (i in 1:t) {
    if (i == 1) {
      comp2_var <- component_var[2]
      w2_i <- 1 / (comp2_var + lfc_var)
      comp2_mean_min <- stats::qnorm(1 - overlap, log2(threshold), sqrt(comp2_var))
      comp2_mu <- max(comp2_mean_min, sum(comp2_pd * w2_i * x) / sum(comp2_pd * w2_i))

      comp3_var <- component_var[3]
      w3_i <- 1 / (comp3_var + lfc_var)
      comp3_mean_max <- stats::qnorm(overlap, -log2(threshold), sqrt(comp3_var))
      comp3_mu <- min(comp3_mean_max, sum(comp3_pd * w3_i * x) / sum(comp3_pd * w3_i))

      comp1_var <- comp1_var_max

      comp1_mu <- 0
    } else {
      comp2_var <- max(0.005, sum(comp2_pd * ((w2_i)^2) * (((x - comp2_mu)^2) - lfc_var)) / (sum(comp2_pd * w2_i^2)))
      w2_i <- 1 / (comp2_var + lfc_var)
      comp2_mean_min <- stats::qnorm((1 - overlap), log2(threshold), sqrt(comp2_var))
      comp2_mu <- max(comp2_mean_min, sum(comp2_pd * w2_i * x) / sum(comp2_pd * w2_i))

      comp3_var <- max(0.005, sum(comp3_pd * ((w3_i)^2) * (((x - comp3_mu)^2) - lfc_var)) / (sum(comp3_pd * w3_i^2)))
      w3_i <- 1 / (comp3_var + lfc_var)
      comp3_mean_max <- stats::qnorm(overlap, -log2(threshold), sqrt(comp3_var))
      comp3_mu <- min(comp3_mean_max, sum(comp3_pd * w3_i * x) / sum(comp3_pd * w3_i))

      comp1_var <- comp1_var_max
      comp1_mu <- 0
    }
  }


  comp1_alpha <- max(comp1_n / sum(set), lower_bound)
  comp2_alpha <- max(comp2_n / sum(set), lower_bound)
  comp3_alpha <- max(comp3_n / sum(set), lower_bound)
  comp4_alpha <- max(comp4_n / (length(set) - sum(set)), lower_bound)
  comp5_alpha <- max(comp5_n / (length(set) - sum(set)), lower_bound)
  comp6_alpha <- max(comp6_n / (length(set) - sum(set)), lower_bound)

  list(
    "mu" = c(comp1_mu, comp2_mu, comp3_mu),
    "var" = c(comp1_var, comp2_var, comp3_var),
    "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha, comp6_alpha)
  )
}


lower_bound <- 0.000000005


#'#' This is data to be included in my package
#'
#' @name data-test.set
#' @docType data
#' @description
#' A list containing a dataframe of genes, effect sizes and effect size variances, and a list with a single gene set.
#' Access using test.set$postdata and test.set$gs respectively. This data is used to carry out tests of the mrema function.
#'
#' @keywords data
"test.set"


#'#' This is data to be included in my package
#'
#' @name data-em.tests
#' @docType data
#' @description
#' A list containing a dataframe of genes, effect sizes and effect size variances, and a list with a single gene set.
#' Access using em.tests$postdata and em.tests$gs respectively. This data is used to carry out tests of the mrema function.
#'
#' @keywords data
"em.tests"
