# fireSense_IgnitionFit 1.0.2

## Breaking changes

- Removed the non-xgboost modelling path entirely (`runGLMMAdaptiveWithSimplifications()`,
  `runGLM.NB()`, `checkData()`, `messageFormulaFn()`, and the formula/AIC code in `buildModel()`),
  together with the disabled time-ordered branch of `runXGBOOST()`.
  `buildModel()` now stops unless `modelAlgorithm` contains `"xgb"`.
- Removed defunct parameters: `escapeFamily`, `ignitionFamily`, `plot_fuelBiomassPerPrediction`,
  `.studyAreaName`, `.plotInitialTime`, `.saveInitialTime`, `.saveInterval`.
- Removed defunct inputs: `fireSense_ignitionFormula`, `climateVariablesForFire`.
- `modelList` no longer carries a `family` element: it was always the `stats::family` function
  itself, never a family object, and nothing read it.

## Bug fixes

- The digest of the covariates is now actually passed to `runXGBOOST()`. The module asked
  `prepareCovariatesOuter()` for `digestOfData$fireSense_ignitionCovariates`, which is not an
  element of what that function returns, so `NULL` was passed and the data played no part in the
  cache key of the per-fold models. **This changes the cache keys: existing caches of the fit are
  invalidated and the models will be refitted once.**

# fireSense_IgnitionFit 1.0.1

First release from `development` since `master` was last updated (2020-04-02). Full history: https://github.com/PredictiveEcology/fireSense_IgnitionFit/compare/8ea03bc...v1.0.1

## Breaking changes

- Removed input `dataFireSense_IgnitionFit` (data.frame).
- Removed parameters: `cores`, `data`, `family`, `formula`, `iterDEoptim`, `iterNlminb`, `lb`, `nlminb.control`, `plot`, `start`, `trace`, `ub`.

## New features

- New inputs: `climateVariablesForFire`, `fireSense_ignitionCovariates`, `fireSense_ignitionFormula`, `ignitionFitRTM`.
- New output: `fireSense_EscapeFitted`.
- New parameters: `.plotInitialTime`, `.plots`, `.seed`, `.studyAreaName`, `crossValType`, `escapeFamily`, `ignitionFamily`, `modelAlgorithm`, `plot_fuelBiomassPerPrediction`, `rescaleVars`, `whichProcessesToFit`.

## Dependencies

- No longer depends on `DEoptim`.
- Now depends on `RhpcBLASctl`, `SHAPforxgboost`, `SpaDES.core`, `caret`, `data.table`, `fireSenseUtils`, `ggplot2`, `ggpubr`, `glmmTMB`, `lightgbm`, `pROC`, `parallelly`, `pemisc`, `reproducible`, `terra`, `xgboost`.

## Testing

- testthat suite and CI (`testthat-module`), including a snapshot of the module's inputs, outputs and parameters in `tests/testthat/test-metadata.R`.
