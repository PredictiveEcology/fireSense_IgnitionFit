---
title: "fireSense_IgnitionFit Manual"
subtitle: "v.1.0.2.9000"
date: "Last updated: 2026-09-24"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
bibliography: citations/references_fireSense_IgnitionFit.bib
link-citations: true
always_allow_html: true
pkgdown:
  as_is: true
---

# fireSense_IgnitionFit Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:fireSense-IgnitionFit) *fireSense_IgnitionFit*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Eliot McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut, cre], Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut], Jean Marchal <jean.d.marchal@gmail.com> [aut], Alex M Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

Fits models of fire ignition and, optionally, fire escape from climate and fuel covariates.
The fitted models are used by *fireSense_IgnitionPredict* and *fireSense_dataPrepPredict*.
The module does not prepare data; *fireSense_dataPrepFit* supplies its inputs.
Earlier versions fitted the piecewise regression models of @Marchal:2017a, @Marchal:2017b and @Marchal:2019.

### What it does

For each process in `whichProcessesToFit` (`"ignition"`, `"escape"`):

1. Covariates are centred and scaled (`fireSenseUtils::prepareCovariatesOuter()`), if `rescaleVars = TRUE`.
   Ignition uses `fireSense_ignitionCovariates`; escape uses `fireSense_escapeCovariates` and only its rows with at least one ignition.
2. Rows are split into 5 cross-validation folds, stratified on whether the response is positive.
   The fit stops if any fold would train without a positive observation.
3. One `xgboost` model (Tweedie objective) is fitted per fold, with the fold as the evaluation set.
   The AUC of each fold is calculated on its held-out rows and the mean is printed.
4. The per-fold models are returned together with the values needed to predict from them (see outputs).

Only `modelAlgorithm = "xgboost"` works: the non-xgboost modelling path was removed, and `buildModel()` now stops for any other algorithm.
`crossValType` does not change the fit; it is only used in the plot filename.
Fits are cached with `reproducible::Cache()`.

### Events

- `init`: schedules `checkData` and `run` at `.runInitialTime`.
- `checkData`: stops if `ignitionFitRTM` lacks the `nonNAs` attribute, or if `whichProcessesToFit` names neither process.
- `run`: fits and, if `.plots` is set, plots. Repeats every `.runInterval` if that is not `NA`.

### Module inputs and parameters

Table \@ref(tab:moduleInputs-fireSense-IgnitionFit) lists the declared inputs.
Fitting escape also needs `fireSense_escapeCovariates` (from *fireSense_dataPrepFit*), which is not declared in the metadata.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-fireSense-IgnitionFit)(\#tab:moduleInputs-fireSense-IgnitionFit)List of (ref:fireSense-IgnitionFit) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
   <th style="text-align:left;"> sourceURL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> fireSense_ignitionCovariates </td>
   <td style="text-align:left;"> data.frame </td>
   <td style="text-align:left;"> Table of aggregated ignition covariates with annual `ignitions` counts, one row per `pixelID` and `year`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ignitionFitRTM </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Template raster at the resolution and extent of `fireSense_ignitionCovariates`. Its resolution and its number of non-NA cells, which must be in the attribute `nonNAs` (`attributes(ignitionFitRTM)$nonNAs`), are stored in the fitted objects. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Parameters are in Table \@ref(tab:moduleParams-fireSense-IgnitionFit).

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-fireSense-IgnitionFit)(\#tab:moduleParams-fireSense-IgnitionFit)List of (ref:fireSense-IgnitionFit) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramClass </th>
   <th style="text-align:left;"> default </th>
   <th style="text-align:left;"> min </th>
   <th style="text-align:left;"> max </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> crossValType </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> time-ord.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Has no effect on the fit: `runXGBOOST()` always uses k-fold cross validation. The first element is only used in the plot filename. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rescaleVars </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> If `TRUE`, covariates are centred and scaled with `scale()` before fitting. The centring and scaling values are returned in the `scaleData` element of the outputs. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> whichProcessesToFit </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> ignition.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> which processes to fit: ignition, escape, or both (the default) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> modelAlgorithm </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> xgboost </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Model type. Only `xgboost` (any value containing 'xgb') works. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plots </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> screen </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> See `?Plots`. If set, plots the predicted response against each covariate, for climate and fuel covariates separately, and saves it as png in `figurePath(sim)`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .runInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> when to start this module? By default, the start time of the simulation. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .runInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. Interval between two runs of this module, expressed in units of simulation time. By default, NA, which means that this module only runs once per simulation. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .seed </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Named list of seeds to use for each event (names). E.g., `list('init' = 123)` will `set.seed(123)` at the start of the init event and unset it at the end. Defaults to `NULL`, meaning that no seeds will be set. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .studyAreaName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Human-readable name for the study area used. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant. </td>
  </tr>
</tbody>
</table>

### Plotting

If `.plots` is set, one figure per process: the predicted response against each covariate, with all other covariates at their mean, for climate and fuel covariates in separate panels.
It is saved as png in `figurePath(sim)`.

### Saving

Nothing is saved apart from the figure.

### Module outputs

Outputs are in Table \@ref(tab:moduleOutputs-fireSense-IgnitionFit).
`fireSense_IgnitionFitted$modelList$model` is the list of per-fold models (`Fold1`, ...) plus `rocs`, the ROC curve of each fold.
`lambdaRescaleFactor` is the number of rows in the covariates divided by the number of non-NA cells in `ignitionFitRTM`.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-fireSense-IgnitionFit)(\#tab:moduleOutputs-fireSense-IgnitionFit)List of (ref:fireSense-IgnitionFit) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> fireSense_EscapeFitted </td>
   <td style="text-align:left;"> fireSense_EscapeFit </td>
   <td style="text-align:left;"> List of `modelList` and `scaleData`, as `fireSense_IgnitionFitted`, with `modelList` of class `fireSense_EscapeFit`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_IgnitionFitted </td>
   <td style="text-align:left;"> fireSense_IgnitionFit </td>
   <td style="text-align:left;"> List of `modelList` (class `fireSense_IgnitionFit`: `model`, the per-fold xgboost models and their ROC curves, `fittingRes`, `lambdaRescaleFactor`, `rescales`) and `scaleData` (centre and scale used to standardise the covariates). </td>
  </tr>
</tbody>
</table>

### Usage


``` r
## in the same `simInit()` call as fireSense_dataPrepFit, which creates the inputs
modules <- c("fireSense_dataPrepFit", "fireSense_IgnitionFit")
params <- list(
  fireSense_IgnitionFit = list(whichProcessesToFit = c("ignition", "escape"), .plots = "png")
)

## after `spades()`
sim$fireSense_IgnitionFitted$modelList$model$Fold1
```

### Links to other modules

- [fireSense_dataPrepFit](https://github.com/PredictiveEcology/fireSense_dataPrepFit) creates the inputs.
- [fireSense_dataPrepPredict](https://github.com/PredictiveEcology/fireSense_dataPrepPredict) and [fireSense_IgnitionPredict](https://github.com/PredictiveEcology/fireSense_IgnitionPredict) use the outputs.

### Getting help

- <https://github.com/PredictiveEcology/fireSense_IgnitionFit/issues>

## References

<!-- autogenerated from bibligraphy -->

