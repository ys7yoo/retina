# Retina Project (collaboration with CBNU College of Medicine)

## Core functions
* [`load_data.m`](load_data.m) - load random stimulus and spike time info from mat files (and save them to txt files)
* [`calc_STA_and_STC.m`](calc_STA_and_STC.m) - calculate STA and STC
* [`find_significant_eigen_values.m`](find_significant_eigen_values.m) - find significant eigen values of STC using a nested hypothesis test

## How to run STA and STC
* See [`demo_sta_vs_stc.m`](demo_sta_vs_stc.m)

## [FUTURE WORK] To add GLM analysis
1. Install [GLMspiketools](https://github.com/ys7yoo/GLMspiketools) developed by [J. Pillow](https://github.com/pillowlab/GLMspiketools).
2. Run `estimGLM.m` to calc STA, GLM, and bilinear GLM for some chosen channels.
