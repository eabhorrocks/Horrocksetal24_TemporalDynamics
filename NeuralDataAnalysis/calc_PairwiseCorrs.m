function rsc = calc_PairwiseCorrs(spikes1, spikes2)

% Uses the method of Kohn and Smith, 2005, simple
% pearson's correlation between two neurons responses to repeats of the 
% same stimulus. Sig. faster than corr(x,y).

% force spike counts to column vectors
spikes1 = spikes1(:); spikes2 = spikes2(:);

% calculate r_sc
rsc = (mean(spikes1.*spikes2) - mean(spikes1)*mean(spikes2))...
    /(std(spikes1).*std(spikes2));



end