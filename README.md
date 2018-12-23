This repository contains analysis code used in Pastoll et al.

We recommend using RStudio (https://www.rstudio.com/) for running and viewing code.

**Suggested set up:**
First run the files set-up.Rmd, Functions.Rmd and LoadData.Rmd. Install any packages that are required for these files to run correctly.

Analyses for each figure are in their own .rmd files. Once libraries are activated, functions initialised and data loaded, then the analysis files can be run independently. Within each .rmd file code should be executed in order.

**Generate a html document to view the code:**
The code can also be viewed as a bookdown document, with each chapter containing analysis for the corresponding figure in the manuscript. This is made possible using the **bookdown** package (https://github.com/rstudio/bookdown). Install the package from from CRAN or Github:

```{r eval=FALSE}
install.packages("bookdown")
# or the development version
# devtools::install_github("rstudio/bookdown")
```
To compile the book then use:
```
bookdown:::serve_book()
```

Please see the page "Get Started" at https://bookdown.org/ for more information.
