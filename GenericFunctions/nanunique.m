function [u, ia, ic] = nanunique(varargin)
x = varargin{1};
t = rand;
while any(x(:)==t)
    t = rand;
end
x(isnan(x)) = t;
[u, ia, ic] = unique(x,varargin{2:end});
u(u==t)=nan;
end