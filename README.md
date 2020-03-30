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

### Web Application:

A (minimal) `shiny` web application to use CIMICE is available [here](https://redsnic.shinyapps.io/CIMICE_WEB/).