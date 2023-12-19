function [p,H,p2] = schreirer_hare(data,grp)

[R,tie] = tiedrank(data);

% compute the anova SS's on the ranks:
[p1,p2] = anovan(R,grp,'display','off','model','interaction');
nvar = size(p2,1)-2;
for l=1:nvar
df(l) = p2{l+1,3};
SS(l) = p2{l+1,2};
end
% total MS
MS = sum(SS)/sum(df);
% compute ratios
for l=1:nvar-1
H(l) = SS(l)/MS;
end

% these are compared to a chi-2-distribution with df deg of freedom:
for l=1:nvar-1
p(l) = 1-chi2cdf(H(l),df(l));
end