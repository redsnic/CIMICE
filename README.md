# CIMICE-R

CIMICE is a tool in the field of tumor phylogenetics and
its goal is to build a Markov Chain (called Cancer Progression Markov Chain, CPMC) in order to model tumor subtypes evolution.
The input of CIMICE is a Mutational Matrix, so a boolean matrix representing altered genes in a
collection of samples. These samples are assumed to be obtained with single-cell DNA analysis techniques and
the tool is specifically written to use the peculiarities of this data for the CMPC construction.
See [this repository](https://github.com/redsnic/tumorEvolutionWithMarkovChains/tree/master/GenotypeEvolutionPaths) for the original Java version of this tool.

### Installation from GitHub:

```{R}
# install package
devtools::install_github("redsnic/CIMICE", build_vignettes = TRUE, dependencies = TRUE)
library(CIMICE)

# simple example
quick_run(example_dataset())
# show guide
browseVignettes("CIMICE")
```

remember to install [**pandoc**](https://pandoc.org/installing.html) to build the vignettes. 

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
