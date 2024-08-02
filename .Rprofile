##| configure reticulate --------------------------------------------------------
##|

# library(reticulate)
# reticulate::use_python("/opt/anaconda3/bin/python")
# reticulate::py_config()

if (benchmarkme::get_cpu()$model_name == "Apple M1 Max") {
  Sys.setenv(RETICULATE_PYTHON = "/Users/sphache/miniconda3/envs/rstudio/bin/python")
} else {
  Sys.setenv(RETICULATE_PYTHON = "/opt/anaconda3/bin/python")
}

library(reticulate)
reticulate::py_config()


##| load libraries --------------------------------------------------------
##|

library(magrittr)
library(usethis)
library(fsnulls)
library(ggplot2)

##| load all functions --------------------------------------------------------
##|

library(devtools)
devtools::load_all()
