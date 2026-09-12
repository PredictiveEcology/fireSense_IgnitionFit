---
title: "fireSense_IgnitionFit Manual"
subtitle: "v.1.0.1.9000"
date: "Last updated: 2026-09-12"
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

This module fits a statistical model to estimate the contributions of climate and fuel to fire ignition.
<!-- TODO -->
Estimate fire ignition (TODO: fill this out) - [@Marchal:2017a; @Marchal:2017b; @Marchal:2019]

### Module summary

Fit a model of fire ignition from climate and fuel covariates.

### Module inputs and parameters

Describe input data required by the module and how to obtain it (e.g., directly from online sources or supplied by other modules) .

Table \@ref(tab:moduleInputs-fireSense-IgnitionFit) shows the full list of module inputs.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-fireSense-IgnitionFit)(\#tab:moduleInputs-fireSense-IgnitionFit)List of (ref:fireSense_IgnitionFit) input objects and their description.</caption>
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
   <td style="text-align:left;"> climateVariablesForFire </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> The column name in the `fireSense_ignitionCovariates` that is climate, in a named list, .e.g. `climateVariablesForFire = list('ignition' = 'MDC')` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_ignitionCovariates </td>
   <td style="text-align:left;"> data.frame </td>
   <td style="text-align:left;"> table of aggregated ignition covariates with annual ignitions </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ignitionFitRTM </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A (template) raster with information with regards to the spatial resolution and geographical extent of `fireSense_ignitionCovariates`. Used to pass this information onto `fireSense_ignitionFitted` Needs to have number of non-NA cells as attribute: (`ignitionFitRTM@data@attributes$nonNAs`), and optionally, `ignitionFitRTM@data@attributes$meanForestB` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_ignitionFormula </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> formula - as a character - describing the model to be fitted. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Summary of user-visible parameters (Table \@ref(tab:moduleParams-fireSense-IgnitionFit))

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
   <td style="text-align:left;"> How the cross validation should happen, time-ordered or regular k-fold crossValidation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> escapeFamily </td>
   <td style="text-align:left;"> function.... </td>
   <td style="text-align:left;"> binomial.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> a family function (must be wrapped with `quote()`) or a character string naming a family function. Only the negative binomial has been implemented For additional details see `?family`. This was formerly `quote(MASS::negative.binomial(theta = 1, link = 'identity'))`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ignitionFamily </td>
   <td style="text-align:left;"> function.... </td>
   <td style="text-align:left;"> poisson, log </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> a family function (must be wrapped with `quote()`) or a character string naming a family function. Only the negative binomial has been implemented For additional details see `?family`. This was formerly `quote(MASS::negative.binomial(theta = 1, link = 'identity'))`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plot_fuelBiomassPerPrediction </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> when generating plots of climate x fuel class, the log of biomass (g/m2) for which to generate predictions across a gradient of climate values. If supplied, it will override any values in `sim$ignitionFitRTM$meanForestB` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rescaleVars </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Attempt to rescale variables? If `rescalers` is defined, use it to rescale variables as `var / rescalers['var']`. Otherwise, `scale()` will be used to rescale variables to `[0,1]`, if they are not already within this range. </td>
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
   <td style="text-align:left;"> Can be `xgboost`, `glmmtmb`, `glm.nb`, `glmmadaptive`, `glm`; only `xgboost` is supported currently </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plots </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> screen </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> See ?Plots. There are a few plots that are made within this module, if set. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> when to do plot </td>
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
   <td style="text-align:left;"> .saveInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. When to start saving output to a file. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. Interval between save events. </td>
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
   <td style="text-align:left;"> Human-readable name for the study area used. If NA, a hash of `studyAreaLarge` will be used. </td>
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

### Events

<!-- TODO -->
Describe what happens for each event type.

### Plotting

<!-- TODO -->
Write what is plotted.

### Saving

<!-- TODO -->
Write what is saved.

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-fireSense-IgnitionFit)).

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
   <td style="text-align:left;"> A fitted model object of class `fireSense_EscapeFit` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_IgnitionFitted </td>
   <td style="text-align:left;"> fireSense_IgnitionFit </td>
   <td style="text-align:left;"> A fitted model object of class `fireSense_IgnitionFit`. </td>
  </tr>
</tbody>
</table>

### Links to other modules

<!-- TODO: add links to other fireSense modules -->
This model can be used to parameterize the fire ignition component of landscape fire models such as fireSense.

### Getting help

- <https://github.com/PredictiveEcology/fireSense_IgnitionFit/issues>

## References

<!-- autogenerated from bibligraphy -->

