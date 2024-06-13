# Flexible neural population dynamics govern the speed and stability of sensory encoding in mouse visual cortex

Please cite the paper if you use the code in this repo!

# Code structure

- 'processSession_XXX' processes the 'basic.mat' file for each session.
- These processed files are then loaded to perform further analysis, plots and compute statistics using the relevant scripts (see each folder)

Note: some of the processSession functions may take over an hour to run. We encourage you to use a parfor loop in the batchProcessSession script to analyse sessions in parallel.

# Instructions
Simply download this repository and the zipped data file (XXX_basic.mat files for each session).

# Requirements
Tested using Matlab 2022a
Required toolboxes (customised copies included)
- DataHigh (Factor analysis)
- CellExplorer (ACG analysis)
- Umakantha et al., 2021 (Population metrics)

