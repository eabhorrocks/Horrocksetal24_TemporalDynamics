function [H, overlap] = overlapCI(means, CIs)


nVals = numel(means);

if nVals ~= numel(CIs)
    error('must provide equal number of means and CIs')
end

overlap = nan(nVals);
H = overlap;

for ival1 = 1:numel(means)
    for ival2 = ival1+1:numel(means)
        
        temp_means = means([ival1,ival2]);
        temp_CIs = CIs([ival1,ival2]);
        
        [~, idx_max] = max(temp_means);
        [~, idx_min] = min(temp_means);
        
        overlap(ival1,ival2) = (temp_means(idx_max)-temp_CIs(idx_min));
        
        H(ival1, ival2) = overlap(ival1,ival2) >= 0;
        
    end
end