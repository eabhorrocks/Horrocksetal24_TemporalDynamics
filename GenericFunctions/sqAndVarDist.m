function [dist, sqDist, varDist] = sqAndVarDist(x,y)

if size(x)~=size(y)
    error('x and y must be equal size')
end

distVec = x-y;
sqDist = sum(distVec.^2);
varDist = numel(x)*var(distVec);
dist = sqDist + varDist;

end