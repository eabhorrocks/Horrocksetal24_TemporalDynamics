function [proportion, SD] = prop(x, includeNaN)
% convert array x to a logical and then calculate proportion of logical
% true values.
if nargin <2 || ~includeNaN % this doesn't work at the moment for == ...
x = x(find(~isnan(x)));
end
x = logical(x);



proportion = sum(x(:))/numel(x);


SD = sqrt((proportion.*(1-proportion))./numel(x));
