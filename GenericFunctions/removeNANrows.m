function X = removeNANrows(X,type)

if nargin>1 && strcmp(type, 'all')
    X(all(isnan(X), 2), :) = [];
else
    X(any(isnan(X), 2), :) = [];
end

end