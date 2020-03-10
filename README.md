This repository contains analysis code used in Pastoll et al., 2020. Inter- and intra-animal variation in the integrative properties of stellate cells in the medial entorhinal cortex, eLife, https://doi.org/10.7554/eLife.52258

We recommend using RStudio (https://www.rstudio.com/) for running and viewing code.

**Some key files:**  
set-up.Rmd  - loads required packages. These packages need to be installed. Additional required packages are listed below.  
Functions.Rmd - initialises functions.  
LoadData.Rmd - loads data.  
ConceptualFigures.Rmd - Figure 1.  
Props_Example.Rmd - Figure 2.  
Props_All.Rmd - Figure 3.  
Props_mixed_model.Rmd - Figure 4.  
Props_FixedEffects.Rmd - Figure 5.  
PCA_mixed_model.Rmd - Figure 7.  

**Suggested set up:**
Make a new RStudio project. Add the files to the project folder. On beggining a new analysis session, first run the files set-up.Rmd, Functions.Rmd and LoadData.Rmd.

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

See the page "Get Started" at https://bookdown.org/ for more information.


Additional packages that need to be installed:  
gridExtra  
kableExtra  
car  
lmerTest
merTools

