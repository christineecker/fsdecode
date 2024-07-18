
#' fs_decode_gls
#'
#' Performs gene expression decoding by using `gls` accounting for spatial autocorrelations under considerations of the different donors.
#'
#' @param fs.overlay FreeSurfer overlay concatenated across hemispheres with NAs in medial wall label(s).
#' @param max.distance DESCRIPTION.
#' @param mirror Either "TRUE" of "FALSE". If true, left hemisphere samples are projected onto the right hemisphere, and vice versa.
#' @param genes Gene list in 'Symbols' to be assessed. If not specified, all abagen genes are used.
#'
#' @export
#'
#' @return data.table with beta, t, p-values, and fdr adjusted p-values for decoded genes.
#'
#' @examples
#' fs.overlay <- c(fsdata::pet.5HT$`5HT1A`$lh, fsdata::pet.5HT$`5HT1A`$rh)
#' fs.overlay[fsdata::fsavg6$medial.wall.verts$both] <- NA
#' genes.decoded <- fs_decode_gls(fs.overlay, max.distance = 5, mirror = F, genes = fsdata::abagen.genes$symbol[1:100])
#'
fs_decode_gls <-
  function(fs.overlay,
           max.distance = 5,
           mirror = FALSE,
           genes = abagen.genes$symbol,
           mc.cores = parallel::detectCores() - 5) {

    if (length(fs.overlay) == 2 * 40962) {
      writeLines("\nUsing fsaverage6 ...")
      fsavg <- "fsavg6"
      mRNA.data <- abagen.mRNA$fsavg6
    } else {
      writeLines("\nUsing fsaverage7 ...")
      fsavg <- "fsavg7"
      mRNA.data <- abagen.mRNA$fsavg7
    }

    ######################### Extract vertex neighborhood ##########################

    sample.neighborhood <- list()
    sample.neighborhood <- get_ahba_sample_vertex_neighborhood(fsaverage = fsavg,
                                                               max.distance = max.distance,
                                                               mirror = mirror)



    ############# Average fs.overlay signal within vertex neighborhood #############

    fs.overlay.sample.data <- sapply(sample.neighborhood, function(x) {mean(fs.overlay[x], na.rm = T)})

    ################################## Fit model ###################################

    writeLines("\nPerforming decoding using gls accounting for spatial autocorrelation ...")

    decoding.results <- list()
    decoding.results <- pbmcapply::pbmclapply(genes, function(x, fs.data) {

      df.gene <- data.frame(
        fs.data,
        mRNA.data[, ..x], # extract by gene name
        mRNA.data$mni_x,
        mRNA.data$mni_y,
        mRNA.data$mni_z,
        mRNA.data$ahba_donor
      )
      names(df.gene) <- c("fs.data", "mRNA", "x", "y", "z", "donor")

      fit <- nlme::gls(
        scale(fs.data) ~ scale(mRNA),
        data = df.gene,
        correlation = nlme::corGaus(form = ~ x + y + z | donor, nugget = T),
        method = "REML"
      )
      #plot(nlme:::Variogram(fit, form = ~ x + y + z | donor, resType = "normalized"))
      coef(summary(fit))[2, c(1, 3, 4)]  ## beta, t-value, p
    },
    mc.cores = mc.cores,
    fs.data = fs.overlay.sample.data) %>%
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
# res <- fs_decode_gls(fs.overlay, max.distance = 5, mirror = F, genes = abagen.genes$symbol[1:10])

