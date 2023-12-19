function [e] = groupedErrorBar(y, err)
hold on,
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    e(i) = errorbar(x, y(:,i), err(:,i), '.');
end
