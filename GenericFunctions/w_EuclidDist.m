function dist = w_EuclidDist(x,y, weights)


dist = sqrt(sum(weights.*(y-x).^2));
