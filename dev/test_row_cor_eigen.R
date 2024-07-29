
##| test row_cor_eigen()
##|

source("R/row_cor_eigen.R")

n.spins <- 1000L
n.genes <- 16000L
n.verts <- 81924L

X <- matrix(runif(n.spins*n.verts, 0, 1), n.spins, n.verts)
dim(X)

Y <- matrix(runif(n.genes*n.verts, 0, 1), n.genes, n.verts)
dim(Y)

r.null <- row_cor_eigen(X, Y)
dim(r.null)


##| test on real data --------------------------------------------------------
##|

data.env <- new.env()
data.env <- fs_create_decode_spins_data_env("both")

fs.overlay <- c(pet.5HT$`5HT1A`$lh, pet.5HT$`5HT1A`$rh)
fs.overlay[fsavg6$medial.wall.verts$both] <- 0

n.spins <- 1000L
spins <- fsnulls::fs_create_spins_bloch(fs.overlay, n.spins)

r.null <- matrix()
r.null <- row_cor_eigen(data.env$mRNA, spins)
dim(r.null)

r.obs <- vector()
r.obs <- as.vector(row_cor_eigen(data.env$mRNA, t(fs.overlay))) #| 1 x n.genes

hist(r.obs)

writeLines("\nComputing permutation p-values ...")
p.perm <- vector()
for (i in 1:nrow(data.env$mRNA)) {
  ## compute two-tailed permutation p value
  p.perm[i] = (sum(abs(r.null[i,]) >= abs(r.obs[i])) + 1) / (1 + n.spins)
}

p.fdr.adj <- p.adjust(p.perm, method = "fdr")

