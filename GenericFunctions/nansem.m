function stdError = nansem(x, dim)
% need to count non-nan elements

if ~exist('dim', 'var') 
    %[~, dim] = max(size(x)); % default to largest dimension!
    dim = 1; % default to first dimension to match matlab funcs
end

nNaN = sum(~isnan(x),dim);

stdError = nanstd(x,[],dim)./sqrt(nNaN);

if isempty(stdError)
    stdError = nan;
end

end