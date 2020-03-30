# CIMICE-R

installation:

```{R}
# install package
devtools::install_github("redsnic/CIMICE-R", subdir = "CIMICE", build_vignettes = TRUE)
library(CIMICE)

# simple example
quick.run(example.dataset())
# show guide
browseVignettes("CIMICE")
```
