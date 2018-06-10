# retina
retina project (collaboration with CBNU College of Medicine)

## Updated Codes (May 28, 2018)
* loadData.m - load data in mat files (and save to txt files)
* calcSTAprestim.m - calcSTA using "pre stim"
* runSTAprestim.m - run STA for a channel


## Codes
* loadData.m - load data in mat files (and save to txt files)
* runSTA.m - run STA for a channel
* runSTAbatch.m - run STA for all the channels
* calcSTA.m - core function to calc STA


## How to run STA
1. Setup expriment number and foldernames in `loadData.m`.
2. Run `runSTAbatch.m` to calc STA for all the channels.


## How to run GLM
1. Install [GLMspiketools](https://github.com/ys7yoo/GLMspiketools) developed by J. Pillow.
2. Run `estimGLM.m` to calc STA, GLM, and bilinear GLM for some chosen channels.
