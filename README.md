This repository contains analysis code used in Pastoll et al.

To run analyses the files Functions.Rmd and LoadData.Rmd must first be run.

The code can be viewed as a bookdown document, with each chapter containing analysis for the corresponding figure in the manuscript. This is made possible using the **bookdown** package (https://github.com/rstudio/bookdown). Install the package from from CRAN or Github:

```{r eval=FALSE}
install.packages("bookdown")
# or the development version
# devtools::install_github("rstudio/bookdown")
```
To compile the book then use:
```bookdown:::serve_book()
```

Please see the page "Get Started" at https://bookdown.org/ for more information.
