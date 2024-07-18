
#' fs_decode_sa_null
#'
#' Decodes FreeSurfer overlay using surrogate maps of the mRNA Principal Components.
#'
#' @param fs.overlay Surface overlay on fsaverage6 to be decoded.
#' @param data.env Data environment with objects `mRNA`, `permuted.pcs`, and `gene.allocation`
#' @param hemi Either 'lh', 'rh', or 'both'
#' @param method Either `pearson`, `kendall`, or `spearman`
#' @param use Either `everything`, `all.obs`, `complete.obs`, `na.or.complete`, or `pairwise.complete.obs`
#'
#' @export
#'
#' @return Returns data.table with observed correlations for each gene, p values, and adjusted p values (fdr, maxT)
#'
#' @examples
#' fs.overlay <- c(fsdata::pet.5HT$`5HT1A`$lh, fsdata::pet.5HT$`5HT1A`$rh)
#' fs.overlay[fsdata::fsavg6$medial.wall.verts$both] <- NA
#' res <- fs_decode_sa_null(fs.overlay, genes = abagen.genes$symbol[1:100])
#' res <- fs_decode_sa_null(fs.overlay)
#'
fs_decode_sa_null <- function(fs.overlay,
                              data.env = NULL,
                              hemi = "both",
                              method = "pearson",
                              use = "pairwise.complete.obs",
                              genes = abagen.genes$symbol) {


  ############################ create data enviroment ############################

  if(is.null(data.env)){
    data.env <- new.env()
    data.env <- fs_create_decode_sa_null_data_env("both", genes)
  }

  ######################## compute observed correlations #########################

  writeLines("\nComputing observed correlations between overlay and all genes ...")
  r.obs <- vector()
  r.obs <- pbapply::pbapply(data.env$mRNA, 1, function(y) { cor(y, fs.overlay, use=use, method = method) } )

  ########################## compute null correlations ###########################

  writeLines("\nComputing null correlations between overlay and all genes ...")
  r.null <- list()
  r.null <- pbapply::pblapply(data.env$permuted.pcs, function(X) {
    apply(X, 1, function(y){ cor(y, fs.overlay, use=use, method = method)})
  })
  r.null <- matrix(unlist(r.null), nrow = length(r.null[[1]]), ncol = length(r.null)) %>% t


  ######################### compute permutation p value ##########################

  writeLines("\nComputing permutation p-values ...")
  gene.allocation <- data.env$gene.allocation
  n.genes <- length(gene.allocation)
  n.perm <- dim(r.null)[2]

  p.perm <- vector()
  for (i in 1:n.genes) {
    ## compute two-tailed permutation p value
    p.perm[i] = (sum(abs(r.null[gene.allocation[i],]) >= abs(r.obs[i])) + 1) / (1 + n.perm)
    ## one-tailed permutation p value H1: R > 0
    #p.perm[i] = (sum(r.pcs.perm.matrix[gene.allocation[i],] >= r.obs[i]) + 1) / (1 + n.perm)
  }

  ######################## compute maxT adjusted p values ########################

  writeLines("\nComputing maxT adjusted p-values ...")
  p.adj <- vector()
  maxR <- apply(abs(r.null), 2, max) # get maximum absolute correlation across PCs within permutations
  for (i in 1:n.genes) {
    p.adj[i] = (sum(maxR >= abs(r.obs[i])) + 1) / (1 + n.perm)
  }

  ######################## compute fdr adjusted p-values #########################

  p.fdr.adj <- p.adjust(p.perm, method = "fdr")

  ######################## create data.table with results ########################

  gene.index <- which(abagen.genes$symbol %in% genes)

  genes.decoded <- data.table::data.table("gene.symbol" = genes,
                                          "gene.entrez" = abagen.genes$entrez[gene.index],
                                          "r" = r.obs,
                                          "p.perm" = p.perm,
                                          "p.maxT.adj" = p.adj,
                                          "p.fdr.adj" = p.fdr.adj)

  return(genes.decoded)
}


################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# fs.overlay <- c(pet.5HT$`5HT1A`$lh, pet.5HT$`5HT1A`$rh)
# fs.overlay[fsavg6$medial.wall.verts$both] <- NA
# res <- fs_decode_sa_null(fs.overlay, genes = abagen.genes$symbol[1:100])
# res <- fs_decode_sa_null(fs.overlay)
