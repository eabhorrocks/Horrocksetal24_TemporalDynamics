function [CI_lower, CI_upper] = medianCI(vals, z)

vals = vals(~isnan(vals));
n = numel(vals);

% get lower and upper indexes based on z and number of elements
lowerCIidx = round(n/2 - (z*sqrt(n))/2);
upperCIidx = round(1+n/2 + (z*sqrt(n))/2);

% sort vals into ascending order.
ordered_vals = sort(vals);

% index into sorted vals to get CIs of median
CI_lower = ordered_vals(lowerCIidx);
CI_upper = ordered_vals(upperCIidx);
