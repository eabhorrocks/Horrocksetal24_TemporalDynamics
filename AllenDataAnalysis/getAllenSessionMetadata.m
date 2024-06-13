%% get allen session metadata

%% load csv file
metadata = readtable('sessionInfo.xlsx');

%% get info about which animals contributed to which state

% state = 0 (stationary), 1 (locomotion), nan (both)
stat_idx = find(metadata.state==0 | isnan(metadata.state));
run_idx = find(metadata.state==1 | isnan(metadata.state));

n_stat_male = sum(strcmp(metadata.sex(stat_idx),'M'))
n_stat_female = sum(strcmp(metadata.sex(stat_idx),'F'))

n_run_male = sum(strcmp(metadata.sex(run_idx),'M'))
n_run_female = sum(strcmp(metadata.sex(run_idx),'F'))