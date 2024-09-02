
# fsdecode

<!-- badges: start -->
<!-- badges: end -->

This package enables gene expression decoding using various models. It relies on the python toolboxes [BrainSMASH](https://brainsmash.readthedocs.io/en/latest/) and [neuromaps](https://netneurolab.github.io/neuromaps/index.html), which have to be installed prior to the `fsdecode` package. Also, it relies on the R package `fsnulls`, which can be installed from [GitHub](https://github.com/christineecker/fsnulls).

## Installation

First, install the `fsnulls` package from [GitHub](https://github.com/christineecker/fsnulls).

Then, specify the path to the python interpreter in the `.Rprofile` file, e.g.,

```{r}
if (benchmarkme::get_cpu()$model_name == "Apple M1 Max") {
  Sys.setenv(RETICULATE_PYTHON = "/Users/sphache/miniconda3/bin/python")
} else {
  Sys.setenv(RETICULATE_PYTHON = "/opt/anaconda3/bin/python")
}

library(reticulate)
reticulate::py_config()
```

or

```r
python.path <- "/opt/anaconda3/bin/python"
reticulate::use_python(python.path)
```

## Downloading the repository

Given the large amount of metadata that is contained within the package, it uses `git annex` to store the data on a remove repository other than github. You therefore also need to install [datalad](https://www.datalad.org), which can be used to clone the repository and to download the data. 

To clone the repository, use

```{bash}
datalad clone https://github.com/christineecker/fsdecode.git
```

To download the data, use

```{bash}
datalad get data/
datalad get inst/extdata
# etc.
```
