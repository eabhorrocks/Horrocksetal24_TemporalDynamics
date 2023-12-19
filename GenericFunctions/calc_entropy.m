function entropy = calc_entropy(prob_dist, correction, nObs)

% Calculate the entropy
temp = - prob_dist .* log2(prob_dist);
temp(~isfinite(temp)) = 0;
entropy = sum(temp(:));

if nargin > 1
    switch correction
        case 'MM'
            m = sum(prob_dist~=0);
            entropy = entropy +(m-1)/(2*nObs);
        case 'BUB'
            % paninski estimation method
    end
end
end