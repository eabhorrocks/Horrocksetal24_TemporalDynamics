function [bestMI, bestbin, p, SSI] = calc_MI(spikeCountCellArray, correction, optimiseBins, sigFlag, nMCSamples, nBinLim)

% 1) Choose number of bins (stimulus and response)
% - same as input, optimise to max info, pre-specified N
% 2) Compute MI using this number of bins, optional correction...
% 3) check significance using shuffle method

% I(X,Y) = H(X) + H(Y) ? H(X,Y).


p = NaN;
allSCs = vertcat(spikeCountCellArray{:});
nObs = numel(allSCs);
rangeSC = min(allSCs):max(allSCs); % range of spike counts
nConds = numel(spikeCountCellArray);
maxBins = numel(rangeSC);

if maxBins>nBinLim
    maxBins = nBinLim;
end

%% optimise bins for response
if optimiseBins
    
    % initialise results for each bin
    MI = zeros(1,maxBins);
    
    for nbins = 2:maxBins
        binedges = quantile(allSCs,nbins-1);
        % gen new tempCellArray
        tempCellArray = cell(1,nConds);
        for icond = 1:nConds
            [~, tempCellArray{icond}] = histc(spikeCountCellArray{icond},[-inf;binedges(:);inf]);
        end
        
        % get p(x,y)
        p_sr = nan*ones(nbins,nConds);
        for icond=1:nConds
            for icount = 1:nbins
                p_sr(icount,icond) = sum(tempCellArray{icond}==icount)/nObs;
            end
        end
        
        p_s = sum(p_sr,1);
        p_r = sum(p_sr,2);
        
        H_S = calc_entropy(p_s,correction,nObs);
        H_R = calc_entropy(p_r,correction,nObs);
        H_SR = calc_entropy(p_sr(:),correction,nObs);
        MI(nbins) = H_S + H_R - H_SR;
        
    end
    
    [~, bestbin] = max(MI);
    
else % no binning of p_r
    
    bestbin = nBinLim;
    
end

%% calc MI using best bin, to use for sig testing
nbins = bestbin;
binedges = quantile(allSCs, nbins-1);

tempCellArray = cell(1,7);
for icond = 1:nConds
    [~, tempCellArray{icond}] = histc(spikeCountCellArray{icond},[-inf;binedges(:);inf]);
end

p_sr = nan*ones(nbins,nConds);
for icond=1:nConds
    for icount = 1:nbins
        p_sr(icount,icond) = sum(tempCellArray{icond}==icount)/nObs;
    end
end

p_s = sum(p_sr,1);
p_r = sum(p_sr,2);

%H_S = calc_entropy(p_s,correction,nObs);
H_S = calc_entropy(p_s,'MLE',nObs); % know probability function perfectly.
H_R = calc_entropy(p_r,correction,nObs);
H_SR = calc_entropy(p_sr(:),correction,nObs);
bestMI = H_S + H_R - H_SR;

%% Stimulus-specific information

n_r = round(p_r*nObs);

% SSI(s) = sumOver(r) p(rIs) * [H_R - H_RIs]

p_rIs = NaN([nbins, nConds]);
p_sIr = NaN([nConds, nbins]);
 for icond = 1:nConds
     for ibin = 1:nbins
     p_rIs(ibin,icond) = sum(tempCellArray{icond}==ibin)/numel(tempCellArray{icond});
     
     p_sIr(icond,ibin) = sum(tempCellArray{icond}==ibin)/...
            n_r(ibin);
     end
 end
 
H_SIr = nan*ones(nbins,1);
for ibin = 1:nbins
    H_SIr(ibin) = calc_entropy(p_sIr(:,ibin), correction,nObs);
end
EntropyTerm = H_S - H_SIr;

SSItemp = NaN([nbins, nConds]);
for icond = 1:nConds 
    for ibin = 1:nbins
    SSItemp(ibin,icond) = (p_rIs(ibin,icond).*EntropyTerm(ibin));
    end
end

SSI = nansum(SSItemp,1);

%% shuffle method 2

if sigFlag
    
    MCInfoVals = NaN([1,nMCSamples]);
    
    % get x and y as column vectors
    x = NaN([nObs,1]);
    xidx = 1;
    y = vertcat(tempCellArray{:});
    for icond = 1:numel(tempCellArray)
        x(xidx:(xidx+numel(tempCellArray{icond})-1)) = repelem(icond,numel(tempCellArray{icond}),1);
        xidx = xidx + numel(tempCellArray{icond});
    end

    
    % Perform the Monte Carlo Trials
    for iSample = 1:nMCSamples
        temp_psr = accumarray({x(randperm(nObs)),y},ones(size(x)))./nObs;
        p_s = sum(temp_psr,2);
        p_r = sum(temp_psr,1);
        H_S = calc_entropy(p_s,correction,nObs);
        H_R = calc_entropy(p_r,correction,nObs);
        H_SR = calc_entropy(p_sr(:),correction,nObs);
        MCInfoVals(iSample) = H_S + H_R - H_SR;
    end
    
    % Calculate the p-value
    p = sum(MCInfoVals>bestMI)/nMCSamples;
    
end




