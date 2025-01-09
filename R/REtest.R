#' Run gene set analysis using Random-Effect permutations
#'
#' @param estimates A three column dataframe with columns corresponding to gene name, lfc and lfcSE, respectively
#' @param lfc_thresh The log2-foldchange above which a gene is considered differentially expressed.
#' @param gene_set A list of named gene sets.
#' @param type Which test to conduct, for differential expression, upregulation, downregulation, "DE", "upreg" and "downreg" respectively
#' @param min.size Numeric value for minimum gene set size, defaults to 5
#'
#' @returns A data.frame with gene set name, set size, proportion of distribution above the threshold and p-value.
#'
#' @export

REtest <- function(estimates, lfc_thresh, gene_set, type = "DE", min.size = 5){
  estimates <- estimates[stats::complete.cases(estimates),]
# testing
  if(type == "DE"){
  weight <- stats::pnorm(-lfc_thresh, estimates$lfc, estimates$lfcSE, lower.tail = TRUE) + stats::pnorm(lfc_thresh, estimates$lfc, estimates$lfcSE, lower.tail = FALSE)
  } else if(type == "upreg"){
    weight <- stats::pnorm(lfc_thresh, estimates$lfc, estimates$lfcSE, lower.tail = FALSE)
  } else if(type == "downreg"){
    weight <- stats::pnorm(-lfc_thresh, estimates$lfc, estimates$lfcSE, lower.tail = TRUE)
  }

  # check gene sets
  gene_set <- lapply(gene_set, function(x) x[x %in% estimates$genes])
  set.size <- lengths(gene_set)
  gene_set <- gene_set[set.size >= min.size]
  set.size <- lengths(gene_set)
  set.weights <- unlist(lapply(gene_set, function(x)mean(weight[estimates$genes %in% x])))
  background.weights <- unlist(lapply(gene_set, function(x)mean(weight[!(estimates$genes %in% x)])))
  enrichment <- set.weights/background.weights


  index <- lapply(gene_set, function(x) which(estimates$genes %in% x))
  # weight perms
  m <- replicate(100000, sample(weight, max(set.size)))
  set.size <- as.character(set.size)
  perm.sizes <- lapply(names(table(set.size)), function(x) {
    x <- as.numeric(x)
    c(mean(colMeans(m[1:x,])), stats::sd(colMeans(m[1:x,])))
  })
  names(perm.sizes) <- names(table(set.size))
  # Assume weights from perms are normally distributed
  p.val <- unlist(lapply(c(1:length(gene_set)), function(x){
    #stats::pnorm(set.weights[x], perm.sizes[[set.size[x]]][1], perm.sizes[[set.size[x]]][2], lower.tail = FALSE)
    1 - truncnorm::ptruncnorm(q = set.weights[x], a = 0, mean = perm.sizes[[set.size[x]]][1], sd = perm.sizes[[set.size[x]]][2])
  }))
  perm.p.val <- unlist(lapply(c(1:length(gene_set)), function(x){
    s <- length(gene_set[[x]])
    mean(colMeans(m[1:s,]) >= set.weights[x])
  }))


  results <- tibble::tibble("Gene.Set" = names(gene_set), "Set.Size" = as.numeric(set.size), "Prop.DE" = set.weights, "Prop.DE.background" = background.weights, "Enrichment" = enrichment, "P.value" = p.val, "Perm.p.value" = perm.p.val, "Index" = index, row.names = NULL)

  # res <- cbind(res, index)
  list("results" = results, "estimates" = estimates)
}
