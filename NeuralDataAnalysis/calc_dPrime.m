%% Calculate d' for neural population spike counts between two stimului
% INPUTS
% - Array of spike counts for two stimuli to compare, of the form:
% (nTrial x nCell) with each row corresponding to a
% trial (observation) and each column to an individual cell (variable)
%
% OUTPUTS
% - d' scalar (unsigned)

% d' = u1 - u2 / (1/2 (var1 + var2))^1/2 (with ref to projected data)

function dprime = calc_dPrime(spCounts1, spCounts2, projType)

u1mean = mean(spCounts1);
u2mean = mean(spCounts2);

switch projType     % generate unit vector to project onto
    case 'diffmeans' 
        vecProj = u1mean - u2mean; % vector through means of two stim responses
        wopt = vecProj/norm(vecProj);
        wopt = wopt';
    
    case 'fld' % fishers linear discriminant (optimal projection?)
        C1 = cov(spCounts1);
        C2 = cov(spCounts2);
        covMatrix = C1+C2; % was prev round(C1+C2,7,'decimals');
        wopt = pinv(covMatrix)*(u1mean-u2mean)'; % more effic. than inv(C1+C2)*(u1-u2)
        wopt = wopt/norm(wopt);
        
    otherwise
        error('Choose valid ''projType'': ''difmeans'' or ''fld''');
end


% for each trial, project the data onto projection unit vector
for i = size(spCounts1,1):-1:1
    s1proj(i,:) = (spCounts1(i,:)*wopt)*wopt;
    s2proj(i,:) = (spCounts2(i,:)*wopt)*wopt;
end

% calculate new means of projected data
meanu1p = mean(s1proj);
meanu2p = mean(s2proj);
diffInpMeans = norm(meanu1p-meanu2p);

s1projVar = var(s1proj);
s2projVar = var(s2proj);
s1var = norm(s1projVar);
s2var = norm(s2projVar);


dprime = (diffInpMeans) ./ sqrt(.5*(s1var + s2var));
