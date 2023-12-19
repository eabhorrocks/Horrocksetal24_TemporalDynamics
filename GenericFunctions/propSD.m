function SD = propSD(p,n)

if numel(n)>1
    error('n must be a single count')
end

SD = sqrt((p.*(1-p))./n);

end