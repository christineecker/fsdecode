
#' fs_decode_lme
#'
#'  Performs neurosynth equivalent decoding of FreeSurfer overlay using a mixed effects model with random slope and intercept.
#'
#' @param fs.overlay FreeSurfer overlay concatenated across hemispheres with NAs in medial wall label(s).
#' @param max.distance DESCRIPTION.
#'
#' @export
#'
#' @return data.table with beta, t, p-values, and fdr adjusted p-values for decoded genes.
#'
#' @examples
#' fs.overlay <- c(fsdata::pet.5HT$`5HT1A`$lh, fsdata::pet.5HT$`5HT1A`$rh)
#' fs.overlay[fsavg6$medial.wall.verts$both] <- NA
#' genes.decoded <- fs_decode_lme(fs.overlay, max.distance = 5, mirror = F)
#'
fs_decode_lme <- function(fs.overlay, max.distance = 5, genes = abagen.genes$symbol) {

  ############################### Define variables ###############################

  if (length(fs.overlay) == 2 * 40962) {
    writeLines("\nUsing fsaverage6 ...")
    fsavg <- "fsavg6"
    mRNA.data <- abagen.mRNA$fsavg6
  } else {
    writeLines("\nUsing fsaverage7 ...")
    fsavg <- "fsavg7"
    mRNA.data <- abagen.mRNA$fsavg7
  }

  donor <- mRNA.data$ahba_donor
  gene.index <- which(names(mRNA.data) %in% genes)
  E <- as.matrix(mRNA.data[, ..gene.index])

  ######################### Extract vertex neighborhood ##########################

  sample.neighborhood <- list()
  sample.neighborhood <- get_ahba_sample_vertex_neighborhood(fsaverage = fsavg,
                                                             max.distance = max.distance,
                                                             mirror = F)

  ############# Average fs.overlay signal within vertex neighborhood #############

  fs.overlay.sample.data <- sapply(sample.neighborhood, function(x) {mean(fs.overlay[x], na.rm = T)})

  ################################## Fit model ###################################

  writeLines("\nPerforming decoding using linear mixed effects model ...")

  decoding.results <- list()
  decoding.results <- pbmcapply::pbmclapply(1:length(genes), function(i) {
    # fit model with random slope and intercept
    fit <- afex::lmer(scale(fs.overlay.sample.data) ~ scale(E[, i]) + (1 + scale(E[, i]) | donor))
    decoding.results <- summary(fit)$coefficients[2, c(1, 4, 5)]
  },
  mc.cores = parallel::detectCores() - 2) %>%
    unlist() %>%
    matrix(., nrow = 3, ncol = length(genes)) %>%
    t

  ####################### Adjust for multiple comparisons ########################

  p.adj <- p.adjust(decoding.results[,3], method = "fdr")

  genes.decoded <- data.table::data.table(
    "gene.symbol" = genes,
    "gene.entrez" = abagen.genes$entrez[which(abagen.genes$symbol %in% genes)],
    "beta" = decoding.results[,1],
    "t"= decoding.results[,2],
    "p.value" = decoding.results[,3],
    "p.fdr.adj" = p.adj
  )

  return(genes.decoded)
}


################################################################################
#                                                                              #
#                                   Testing                                    #
#                                                                              #
################################################################################

# fs.overlay <- c(pet.5HT$`5HT1A`$lh, pet.5HT$`5HT1A`$rh)
# fs.overlay[fsavg6$medial.wall.verts$both] <- NA
# genes.decoded <- fs_decode_lme(fs.overlay, max.distance = 5, genes = abagen.genes$symbol[1:10])
# head(genes.decoded)
# genes.decoded[which(genes.decoded$gene.symbol == "HTR1A")]
