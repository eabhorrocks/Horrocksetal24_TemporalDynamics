function zScored = z_score(vals, mu, sigma)
%% function to z-score values according to a known normal distribution.
% example usage for analysing PSTHs: 
% Generate normal dist by z-scoring baseline period and then use this to
% z-score the rest of the PSTH.

zScored = (vals-mu)/sigma;

end