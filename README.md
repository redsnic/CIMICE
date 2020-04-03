# CIMICE-R

### Installation:

```{R}
# install package
devtools::install_github("redsnic/CIMICE-R", subdir = "CIMICE", build_vignettes = TRUE)
library(CIMICE)

# simple example
quick.run(example.dataset())
# show guide
browseVignettes("CIMICE")
```

### Input Formats:

#### CAPRI format

```
s/g    gene_1 gene_2 ... gene_n
sample_1 1 0 ... 0
...
sample_m 1 1 ... 1
```

#### CAPRIpop format

```
s/g    gene_1 gene_2 ... gene_n freq
sample_1 1 0 ... 0 freq_s1
...
sample_m 1 1 ... 1 freq_sm
```

(With this format, feature selection for most mutated samples or genes is disabled in the web application)

### Web Application:

A simple `shiny` web application to use CIMICE is available [here](https://redsnic.shinyapps.io/CIMICE_WEB/).