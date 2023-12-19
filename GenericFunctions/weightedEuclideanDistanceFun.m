function dist = weightedEuclideanDistanceFun(x,Y, weights)

dist = nan([size(Y,1), 1]);

for irow = 1:size(Y,1)
dist(irow) = sqrt(sum(weights.*(Y(irow)-x).^2));
end