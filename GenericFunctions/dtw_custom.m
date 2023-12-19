function [dist,sqDist,varDist,ix,iy] = dtw_custom(x,y,maxSamp) 

[~, ix, iy] = dtw(x,y,maxSamp);

x = x(ix);
y = y(iy);


distVec = x-y;
sqDist = (distVec)*(distVec)'; %sq euclidean
varDist = numel(ix)*var(distVec);

dist = sqDist+varDist;
