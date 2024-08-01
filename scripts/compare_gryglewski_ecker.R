
##| output dir --------------------------------------------------------
##|

library(ggplot2)

out.dir <- "output/Gryglewski_comparison/"
dir.create(out.dir)

################################################################################
#                                                                              #
#                                    HTR1A                                     #
#                                                                              #
################################################################################

##| Plot overlays --------------------------------------------------------
##|

cm <- fs_plot(Gryglewski$HTR1A$lh, "fsaverage6", "lh", c(3.5, 7))
fs_save_plot_hemi(cm, out.dir, "Gryglewski_HTR1A.png", "lh", "log2(mRNA)")

cm <- fs_plot_gene("HTR1A")
fs_save_plot_hemi(cm, out.dir, "Ecker_HTR1A.png", "lh", "mRNA")

df <- data.frame(
  "Gryglewski" = Gryglewski$HTR1A$lh,
  "Ecker" = fs_extract_gene_mRNA("HTR1A")$lh
)
str(df)

##| Compute spatial correlation --------------------------------------------------------
##|

filename <- system.file("extdata/lh.rh.5HT1A.knn30000.sa1000.rda", package = "fsnulls")
load(filename)

r <- fs_spatial_cor_with_nulls(
  fs.x = df$Ecker,
  fs.y = df$Gryglewski,
  fs.x.nulls = sa.5HT1A[,1:fsavg6$nverts],
  method = "pearson"
)
print(r)

##| Create Scatterplot --------------------------------------------------------
##|

label <- bquote(r[pearson] == .(round(r[[1]], 3)) ~ "***")

pl <- ggplot(df, aes(x = Gryglewski, y = Ecker)) %>%
  apply_theme_nature(base.size = 5, font = "Arial") +
  geom_point(
    shape = 18,
    color = "#ad3436",
    alpha = 0.25,
    size = 0.05
  ) +
  geom_smooth(method='lm', size = 0.1, color = "black", se = TRUE) +
  xlab("log2(mRNA) Gryglewski et al.") + ylab("mRNA Ecker et al.") +
  annotate("text", x = 3, y = Inf, label = label, hjust = 0, vjust = +2, size = 2) +
  xlim(3,8) + ylim(0, 0.8)
pl

ggsave(
  filename = paste0(out.dir, "HTR1A_cor.pdf"),
  plot = pl,
  height = 4.5,
  width = 6,
  units = "cm",
  dpi = 300,
  device = cairo_pdf
)


################################################################################
#                                                                              #
#                                    HTR2A                                     #
#                                                                              #
################################################################################

##| Plot overlays --------------------------------------------------------
##|

cm <- fs_plot(Gryglewski$HTR2A$lh, "fsaverage6", "lh", c(7, 9))
fs_save_plot_hemi(cm, out.dir, "Gryglewski_HTR2A.png", "lh", "log2(mRNA)")

cm <- fs_plot_gene("HTR2A", range = c(0.5, 1))
fs_save_plot_hemi(cm, out.dir, "Ecker_HTR2A.png", "lh", "mRNA")

df <- data.frame(
  "Gryglewski" = Gryglewski$HTR2A$lh,
  "Ecker" = fs_extract_gene_mRNA("HTR2A")$lh
)

##| Compute spatial correlation --------------------------------------------------------
##|

filename <- system.file("extdata/lh.rh.5HT2A.knn10000.sa1000.rda", package = "fsnulls")
load(filename)

r <- fs_spatial_cor_with_nulls(
  fs.x = df$Ecker,
  fs.y = df$Gryglewski,
  fs.x.nulls = sa.5HT2A[,1:fsavg6$nverts],
  method = "pearson"
)

label <- bquote(r[pearson] == .(round(r[[1]], 3)) ~ "***")

##| Create Scatterplot --------------------------------------------------------
##|

pl <- ggplot(df, aes(x = Gryglewski, y = Ecker)) %>%
  apply_theme_nature(base.size = 5, font = "Arial") +
  geom_point(
    shape = 18,
    color = "#ad3436",
    alpha = 0.25,
    size = 0.05
  ) +
  geom_smooth(method='lm', size = 0.1, color = "black", se = F) +
  xlab("log2(mRNA) Gryglewski et al.") + ylab("mRNA Ecker et al.") +
  annotate("text", x = 7.5, y = Inf, label = label, hjust = 0, vjust = +2, size = 2) +
  xlim(7.5,9) + ylim(0.7, 1)
print(pl)

ggsave(
  filename = paste0(out.dir, "HTR2A_cor.pdf"),
  plot = pl,
  height = 4.5,
  width = 6,
  units = "cm",
  dpi = 300,
  device = cairo_pdf
)

################################################################################
#                                                                              #
#                                    HTR4                                      #
#                                                                              #
################################################################################

##| Plot overlays --------------------------------------------------------
##|

cm <- fs_plot(Gryglewski$HTR4$lh, "fsaverage6", "lh", c(4, 5))
fs_save_plot_hemi(cm, out.dir, "Gryglewski_HTR4.png", "lh", "log2(mRNA)")

cm <- fs_plot_gene("HTR4", range = c(0.34, 0.64))
fs_save_plot_hemi(cm, out.dir, "Ecker_HTR4.png", "lh", "mRNA")

df <- data.frame(
  "Gryglewski" = Gryglewski$HTR4$lh,
  "Ecker" = fs_extract_gene_mRNA("HTR4")$lh
)

##| Compute spatial correlation --------------------------------------------------------
##|

filename <- system.file("extdata/lh.rh.5HT4.knn20000.sa1000.rda", package = "fsnulls")
load(filename)

r <- fs_spatial_cor_with_nulls(
  fs.x = df$Ecker,
  fs.y = df$Gryglewski,
  fs.x.nulls = sa.5HT4[,1:fsavg6$nverts],
  method = "pearson"
)
print(r)

##| Create Scatterplot --------------------------------------------------------
##|

label <- bquote(r[pearson] == .(round(r[[1]], 3)) ~ "***")

pl <- ggplot(df, aes(x = Gryglewski, y = Ecker)) %>%
  apply_theme_nature(base.size = 5, font = "Arial") +
  geom_point(
    shape = 18,
    color = "#ad3436",
    alpha = 0.25,
    size = 0.05
  ) +
  geom_smooth(method='lm', size = 0.1, color = "black", se = F) +
  xlab("log2(mRNA) Gryglewski et al.") + ylab("mRNA Ecker et al.") +
  annotate("text", x = 4.1, y = Inf, label = label, hjust = 0, vjust = +2, size = 2) +
  xlim(4.1,5.1) + ylim(0.25, 0.75)

ggsave(
  filename = paste0(out.dir, "HTR4_cor.pdf"),
  plot = pl,
  height = 4.5,
  width = 6,
  units = "cm",
  dpi = 300,
  device = cairo_pdf
)
