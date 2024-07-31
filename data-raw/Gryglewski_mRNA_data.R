## code to prepare `Gryglewski_mRNA_data` dataset goes here

## Note. Data downloaded from http://www.meduniwien.ac.at/neuroimaging/mRNA.html

library(freesurferformats)

data.dir <- "inst/extdata/Gryglewski/"

##| HTR1A --------------------------------------------------------
##|

lh.filename <- paste0(data.dir, "3350_HTR1A/3350_mRNA_mirr_lh.mgh")
rh.filename <- paste0(data.dir, "3350_HTR1A/3350_mRNA_mirr_rh.mgh")

lh.HTR1A <- read.fs.morph(lh.filename)
rh.HTR1A <- read.fs.morph(rh.filename)

lh.HTR1A.fsavg6 <- fs_map2fsavg6(lh.HTR1A, "lh")
rh.HTR1A.fsavg6 <- fs_map2fsavg6(rh.HTR1A, "rh")

lh.HTR1A.fsavg6[fsavg6$medial.wall.verts$lh] <- NA
rh.HTR1A.fsavg6[fsavg6$medial.wall.verts$rh] <- NA

##| HTR2A --------------------------------------------------------
##|

lh.filename <- paste0(data.dir, "3356_HTR2A/3356_mRNA_mirr_lh.mgh")
rh.filename <- paste0(data.dir, "3356_HTR2A/3356_mRNA_mirr_rh.mgh")

lh.HTR2A <- read.fs.morph(lh.filename)
rh.HTR2A <- read.fs.morph(rh.filename)

lh.HTR2A.fsavg6 <- fs_map2fsavg6(lh.HTR2A, "lh")
rh.HTR2A.fsavg6 <- fs_map2fsavg6(rh.HTR2A, "rh")

lh.HTR2A.fsavg6[fsavg6$medial.wall.verts$lh] <- NA
rh.HTR2A.fsavg6[fsavg6$medial.wall.verts$rh] <- NA

##| HTR4 --------------------------------------------------------
##|

lh.filename <- paste0(data.dir, "3360_HTR4/3360_mRNA_mirr_lh.mgh")
rh.filename <- paste0(data.dir, "3360_HTR4/3360_mRNA_mirr_rh.mgh")

lh.HTR4 <- read.fs.morph(lh.filename)
rh.HTR4 <- read.fs.morph(rh.filename)

lh.HTR4.fsavg6 <- fs_map2fsavg6(lh.HTR4, "lh")
rh.HTR4.fsavg6 <- fs_map2fsavg6(rh.HTR4, "rh")

lh.HTR4.fsavg6[fsavg6$medial.wall.verts$lh] <- NA
rh.HTR4.fsavg6[fsavg6$medial.wall.verts$rh] <- NA

##| concatenate to list --------------------------------------------------------
##|

Gryglewski <- list(
  "HTR1A" = list("lh" = lh.HTR1A.fsavg6,
                 "rh" = rh.HTR1A.fsavg6),
  "HTR2A" = list("lh" = lh.HTR2A.fsavg6,
                 "rh" = rh.HTR2A.fsavg6),
  "HTR4" = list("lh" = lh.HTR4.fsavg6,
                 "rh" = rh.HTR4.fsavg6)
)

save(Gryglewski, file = "data/Grygleswski.rda")
