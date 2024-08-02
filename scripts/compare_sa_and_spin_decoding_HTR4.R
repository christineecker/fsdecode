
rm(list = ls()) %>% gc()

################################################################################
#                                                                              #
#                                    HTR1A                                     #
#                                                                              #
################################################################################

################################################################################
## Load Data
################################################################################

fs.overlay <- c(pet.5HT$`5HT4`$lh, pet.5HT$`5HT4`$rh)

##| compute spins --------------------------------------------------------
##|

fs.spins <- fsnulls::fs_create_spins_bloch(fs.overlay, 1000L)

##| save a couple of examples --------------------------------------------------------
##|

write.dir <- "output/Nulls/HTR4/"
fs::dir_create(write.dir)

for (i in 1:5) {
  spin <- sample(1:1000, 1)
  plot <- fs_plot(fs.spins[spin, ], hemi = "lh")

  filename <- paste0("HTR4_spin", spin, ".png")
  fs_save_plot_view_angle(
    plot,
    write.dir,
    filename,
    angle.set = "sd_lateral_lh",
    colorbar = "no",
    legend.label = ""
  )
}

##| load surrogates --------------------------------------------------------
##|

# Note. precomputed using `fsnulls` package

filename <- system.file("extdata/lh.rh.5HT4.knn20000.sa1000.rda", package = "fsnulls")
load(filename)
sa.5HT4[fsavg6$medial.wall.verts$both] <- 0

##| save a couple of examples --------------------------------------------------------
##|

for (i in 1:5) {
  sa <- sample(1:1000, 1)
  plot <- fs_plot(sa.5HT4[sa, ], hemi = "lh")

  filename <- paste0("HTR4_sa", sa, ".png")
  fs_save_plot_view_angle(
    plot,
    write.dir,
    filename,
    angle.set = "sd_lateral_lh",
    colorbar = "no",
    legend.label = ""
  )
}

################################################################################
## Decode in various ways
################################################################################

gene.index <- which(abagen.genes$symbol == "HTR4")
data.env <- new.env()
data.env <- fs_create_decode_spins_data_env("both")
genes <- abagen.genes$symbol

##| Using spin models --------------------------------------------------------
##|

genes.spins <- fs_decode_with_nulls(
  fs.overlay,
  fs.nulls = fs.spins,
  data.env = data.env,
  genes = genes
)
print(genes.spins)
print(genes.spins[gene.index,])

# > print(genes.spins[gene.index,])
# gene.symbol gene.entrez         r      p.perm  p.maxT.adj   p.fdr.adj
# <char>       <num>     <num>       <num>       <num>       <num>
#   1:        HTR4        3360 0.8462412 0.000999001 0.001998002 0.001012538

##| using surrogates --------------------------------------------------------
##|

genes.sa <- fs_decode_with_nulls(
  fs.overlay,
  fs.nulls = sa.5HT4,
  data.env = data.env,
  genes = genes
)
print(genes.sa)
print(genes.sa[gene.index,])

# > print(genes.sa[gene.index,])
# gene.symbol gene.entrez         r    p.perm p.maxT.adj p.fdr.adj
# <char>       <num>     <num>     <num>      <num>     <num>
#   1:        HTR4        3360 0.8462412 0.2117882          1         1

##| using surrogate gradients --------------------------------------------------------
##|

##| Note. uses different data.env ...
rm(data.env) %>% gc() # free up some RAM

genes.sa_null <- fs_decode_sa_null(fs.overlay,
                                   data.env = NULL,
                                   hemi = "both",
                                   method = "pearson",
                                   use = "pairwise.complete.obs",
                                   genes = genes)
print(genes.sa_null)
print(genes.sa_null[gene.index,])

# > print(genes.sa_null[gene.index,])
# gene.symbol gene.entrez         r     p.perm p.maxT.adj p.fdr.adj
# <char>       <num>     <num>      <num>      <num>     <num>
#   1:        HTR4        3360 0.3105343 0.03196803  0.4625375 0.1230927


################################################################################
## Compare Models
################################################################################

##| Compute FDR --------------------------------------------------------
##|

n.genes <- length(abagen.genes$symbol)

fpr.spins <- sapply(c(0.05, 0.01, 0.001), function(p) {
  length(which(genes.spins$p.maxT.adj < p)) / n.genes
})
fpr.spins

fpr.sa <- sapply(c(0.05, 0.01, 0.001), function(p) {
  length(which(genes.sa$p.maxT.adj < p)) / n.genes
})
fpr.sa

fpr.sa_null <- sapply(c(0.05, 0.01, 0.001), function(p) {
  length(which(genes.sa_null$p.maxT.adj < p)) / n.genes
})
fpr.sa_null

##| Concatenate data frame --------------------------------------------------------
##|

df.fpr <- data.frame(
  "p" = factor(rep(c("0.05", "0.01", "0.001"), 3), levels = c("0.05", "0.01", "0.001")),
  "method" = factor(c(rep("Alexander-Bloch-2018", 3), rep("Burt-2020", 3), rep("a-null", 3)),
                    levels = c("Alexander-Bloch-2018", "Burt-2020", "a-null")),
  "FPR" = c(fpr.spins, fpr.sa, fpr.sa_null)
)
df.fpr


##| Create figure --------------------------------------------------------
##|

xlabel <- expression("adjusted p-value (p"[adj]*")")
ylabel <- expression("P(p < p"[adj]*")")

p <- ggplot(data = df.fpr, aes(x = p, y = FPR, group = method, color = method, fill = method)) %>%
  apply_theme_nature(base.size = 5, font = "Arial") +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = c('#568ab7', '#ad3436', 'orange'),
                     labels = c("Alexander-Bloch-2020", "Burt-2018", expression(alpha*"-null"))) +
  xlab(xlabel) + ylab(ylabel) +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.position = c(0.8, 0.9),legend.title=element_blank()) +
  theme(legend.text=element_text(size=5))

print(p)

ggsave(
  filename = paste0(write.dir, "HTR4_cor.pdf"),
  plot = p,
  height = 4.5,
  width = 6,
  units = "cm",
  dpi = 300,
  device = cairo_pdf
)


##| Overlap between gene sets --------------------------------------------------------
##|

p.threshold <- 0.05

df <- list(
  "Burt-2020" = genes.sa$gene.symbol[which(genes.sa$p.maxT.adj < p.threshold)],
  "sa_null" = genes.sa_null$gene.symbol[which(genes.sa_null$p.maxT.adj < p.threshold)]
)

myV <- nVennR::plotVenn(df)

cairo_pdf(
  filename = paste0(write.dir, "Venn_0.05.pdf"),
  width = 5,
  height = 4,
  onefile = TRUE
)
nVennR::showSVG(
  myV,
  opacity = 0.2,
  setColors = c('#ad3436', 'orange'),
  labelRegions = F,
  fontScale = 2,
  borderWidth = 3,
  systemShow = F
)
dev.off()
