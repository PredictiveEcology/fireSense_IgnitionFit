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
