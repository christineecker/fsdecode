

write.dir <- "output/mRNA_both_hemis/"
dir.create(write.dir)

##| Load mRNA data --------------------------------------------------------
##|

if (!exists("mRNA")) {
  mRNA <- load_mRNA_expression_maps()
}

##| Select HTR1A --------------------------------------------------------
##|

gene <- "HTR1A"
range <- c(0.2, 0.7)

for (hemi in c("lh", "rh")) {
  index <- get_gene_index(gene)

  if (hemi == "lh") {
    data <- mRNA$lh.rh.mRNA.fsavg6.fwhm5$lh[index, ]
  } else {
    data <- mRNA$lh.rh.mRNA.fsavg6.fwhm5$rh[index, ]
  }

  fs.plot <- fs_plot(data,
                     hemi = hemi,
                     range = range)
  filename <- paste0(gene, "_", hemi, ".png")
  fs_save_plot_hemi(fs.plot,
                    write.dir,
                    filename,
                    angle.set = hemi,
                    legend.label = "mRNA")
}

##| Select HTR2A --------------------------------------------------------
##|

gene <- "HTR2A"
index <- get_gene_index(gene)
range <- c(0.5, 1)

for (hemi in c("lh", "rh")) {
  index <- get_gene_index(gene)

  if (hemi == "lh") {
    data <- mRNA$lh.rh.mRNA.fsavg6.fwhm5$lh[index, ]
  } else {
    data <- mRNA$lh.rh.mRNA.fsavg6.fwhm5$rh[index, ]
  }

  fs.plot <- fs_plot(data,
                     hemi = hemi,
                     range = range)

  filename <- paste0(gene, "_", hemi, ".png")
  fs_save_plot_hemi(fs.plot,
                    write.dir,
                    filename,
                    angle.set = hemi,
                    legend.label = "mRNA")
}


##| Select HTR4 --------------------------------------------------------
##|

gene <- "HTR4"
index <- get_gene_index(gene)
range <- c(0.4, 0.6)

for (hemi in c("lh", "rh")) {
  index <- get_gene_index(gene)

  if (hemi == "lh") {
    data <- mRNA$lh.rh.mRNA.fsavg6.fwhm5$lh[index, ]
  } else {
    data <- mRNA$lh.rh.mRNA.fsavg6.fwhm5$rh[index, ]
  }

  fs.plot <- fs_plot(data,
                     hemi = hemi,
                     range = range)

  filename <- paste0(gene, "_", hemi, ".png")
  fs_save_plot_hemi(fs.plot,
                    write.dir,
                    filename,
                    angle.set = hemi,
                    legend.label = "mRNA")
}
