
##| Note --------------------------------------------------------
##|
##| This script compares the correlation between two genes
##| prior and post spatial correlation
##|

out.dir <- "output/gene_correlations/"
dir.create(out.dir)

mRNA.vertex <- load_mRNA_expression_maps()$lh.rh.mRNA.fsavg6.fwhm5$lh
mRNA.samples <- abagen.mRNA$fsavg6
str(mRNA.samples)

##| Get index for spatially interpolated data
##|

mRNA.index <- grep("GABR", abagen.genes$symbol)
mRNA.vertex.gab <- mRNA.vertex[mRNA.index,]
mRNA.vertex.gab <- t(mRNA.vertex.gab)

##| Get index for sample data
##|

sample.index <- grep("GABR", names(mRNA.samples))
mRNA.samples.gab <- mRNA.samples[,.SD,.SDcols=sample.index]
str(mRNA.samples.gab)

r.mRNA.samples <- cor(mRNA.samples.gab)
r.mRNA.samples <- r.mRNA.samples[lower.tri(r.mRNA.samples, diag = FALSE)]
r.mRNA.vertex <- cor(mRNA.vertex.gab, use = "pairwise.complete.obs")
r.mRNA.vertex <- r.mRNA.vertex[lower.tri(r.mRNA.vertex, diag = FALSE)]

r <- cor(r.mRNA.samples, r.mRNA.vertex, method = "spearman")

plot(r.mRNA.samples, r.mRNA.vertex)


##| Create Scatterplot --------------------------------------------------------
##|

df <- data.frame("Samples" = r.mRNA.samples, "Vertices" = r.mRNA.vertex)
plot(df$Samples, df$Vertices)

label <- bquote(r[spearman] == .(round(r[[1]], 3)) ~ "***")

pl <- ggplot(df, aes(x = Samples, y = Vertices)) %>%
  apply_theme_nature(base.size = 5, font = "Arial") +
  geom_point(
    shape = 19,
    color = "#ad3436",
    alpha = 0.5,
    size = 1
  ) +
  geom_smooth(method='lm', size = 0.3, color = "grey25", se = F) +
  xlab("r(genes) - AHBA resolution") + ylab("r(genes) - vertex resolution") +
  annotate("text", x = -0.5, y = Inf, label = label, hjust = 0, vjust = +2, size = 2)
pl

ggsave(
  filename = paste0(out.dir, "gene_correlations.pdf"),
  plot = pl,
  height = 5,
  width = 7,
  units = "cm",
  dpi = 300,
  device = cairo_pdf
)


##| Create correlation heatmaps --------------------------------------------------------
##| for samples
##|

r.mRNA.samples <- cor(mRNA.samples.gab)

## Heatmap Colorscale
col.map <- colorRampPalette(rev(c("#fabb3e", "white", "darkgrey")))(20)

p <- pheatmap::pheatmap(
  r.mRNA.samples,
  color = col.map,
  border_color = "white",
  cellwidth = 14,
  cellheight = 10,
  cluster_rows = F,
  cluster_cols = F,
  method = "complete",
  display_numbers = T,
  fontsize = 4,
  number_color = "grey40",
  angle_col = 45,
  legend = T,
  breaks = seq(-0.9, 0.9, 0.1)
)%>%
  ggplotify::as.grob()

ggsave(filename = paste0(out.dir, "pheat_samples.pdf"),
       p,
       units = "cm",
       height = 8,
       width = 9)


##| Create correlation heatmaps --------------------------------------------------------
##| for vertices
##|

r.mRNA.vertex <- cor(mRNA.vertex.gab, use = "pairwise.complete.obs")

## Heatmap Colorscale
col.map <- colorRampPalette(rev(c("#fabb3e", "white", "darkgrey")))(20)

p <- pheatmap::pheatmap(
  r.mRNA.vertex,
  color = col.map,
  border_color = "white",
  cellwidth = 14,
  cellheight = 10,
  cluster_rows = F,
  cluster_cols = F,
  method = "complete",
  display_numbers = T,
  fontsize = 4,
  number_color = "grey40",
  angle_col = 45,
  legend = T,
  breaks = seq(-0.9, 0.9, 0.1),
)%>%
  ggplotify::as.grob()

ggsave(filename = paste0(out.dir, "pheat_vertices.pdf"),
       p,
       units = "cm",
       height = 8,
       width = 9)
