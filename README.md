# README

This repository contains an R script and data to reproduce the results presented in "Estimating the Impact of Cycling Infrastructure Improvement on Cycling Infrastructure Usage: A Fixed-Effect Spatial Lag Model Approach."

## How to reproduce the results

You need R version 4.3.2. Recreate the `renv` environment using the `renv.lock` file.

In your R console, execute the following:

```R
renv::restore()
source("script.R", encoding = "UTF-8")
```
