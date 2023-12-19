function stdError = sem(x, dim)

if ~exist('dim', 'var') 
    %[~, dim] = max(size(x)); % default to largest dimension!
    dim = 1; % default to first dim
end

stdError = std(x,[],dim)/sqrt(size(x,dim));

if isempty(stdError)
    stdError = nan;
end

end