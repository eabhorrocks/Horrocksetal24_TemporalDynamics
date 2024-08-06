# Flexible neural population dynamics govern the speed and stability of sensory encoding in mouse visual cortex

If you find this code useful for your research, please cite the paper:
Horrocks, E.A.B., Rodrigues, F.R. & Saleem, A.B. Flexible neural population dynamics govern the speed and stability of sensory encoding in mouse visual cortex. Nat Commun 15, 6415 (2024). https://doi.org/10.1038/s41467-024-50563-y

# Code structure

- 'processSession_XXX' processes the 'basic.mat' file for each session.
- These processed files are then loaded to perform further analysis, plots and compute statistics using the relevant scripts (see each folder)

Note: some of the processSession functions may take over an hour to run. We encourage you to use a parfor loop in the batchProcessSession script to analyse sessions in parallel.

# Instructions
Simply download this repository and the zipped data file (XXX_basic.mat files for each session).
Dataset available at https://figshare.com/articles/dataset/Data_for_Horrocks_et_al_2024/26031226/1

# Requirements
Tested using Matlab 2022a
Required toolboxes (customised copies included)
- DataHigh (Factor analysis)
- CellExplorer (ACG analysis)
- Umakantha et al., 2021 (Population metrics)

