---
title: Setup
---

Please follow the steps below and install the required software **before** the scheduled workshop.

<!--
FIXME: Setup instructions live in this document. Please specify the tools and
the data sets the Learner needs to have installed.

## Data Sets

FIXME: place any data you want learners to use in `episodes/data` and then use
       a relative link ( [data zip file](data/lesson-data.zip) ) to provide a
       link to it, replacing the example.com link.
Download the [data zip file](https://example.com/FIXME) and unzip it to your Desktop
-->
## RStudio Setup

We use RStudio for coding in R.

[Click here and follow the instructions](https://posit.co/download/rstudio-desktop/) to install RStudio Desktop in your system.

::::::::::::::::: discussion

### R packages

Most workshops using R will require the installation of specific packages. Make sure to check in advance with the workshop organisers what packages need to be installed. 

You can install packages from CRAN using:

```r
install.packages("package_name")
```

If your package is in a different R repository, such as Bioconductor or GitHub, you may need the [BiocManager](https://www.bioconductor.org/install/) or [devtools](https://devtools.r-lib.org/) packages to install them. To install BiocManager:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

For devtools, you can simply do:

```r
install.packages("devtools")
```

You can then install packages directly from GitHub with:

```r
devtools::install_github("username/reponame")
```

::::::::::::::::::::::::::::


<!--
READ HERE FOR OS-SPECIFIC INSTRUCTIONS

Setup for different systems can be presented in dropdown menus via a `spoiler`
tag. They will join to this discussion block, so you can give a general overview
of the software used in this lesson here and fill out the individual operating
systems (and potentially add more, e.g. online setup) in the solutions blocks.
-->

<!--
:::::::::::::::: spoiler

### Windows

Use PuTTY

::::::::::::::::::::::::

:::::::::::::::: spoiler

### MacOS

Use Terminal.app

::::::::::::::::::::::::


:::::::::::::::: spoiler

### Linux

Use Terminal

::::::::::::::::::::::::
-->

