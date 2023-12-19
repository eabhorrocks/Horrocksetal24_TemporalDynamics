%% get noise Corrs

function [noiseCorr, signalCorr, totalCorr] = calc_noiseCorr(spikes1, spikes2, nShuffle)

spikes1=spikes1(:);
spikes2=spikes2(:);

if numel(spikes1)~=numel(spikes2)
    error('spikes1 and spikes2 must be vectors of spike counts with equal length')
end

if nargin<3
    nShuffle = numel(spikes1);
end

nCounts = numel(spikes1);


%% calculate total correlation

% totalCorr = (mean(spikes1.*spikes2) - mean(spikes1)*mean(spikes2))...
%     /(std(spikes1).*std(spikes2));

%totalCorr = corr(spikes1,spikes2,'type', 'Pearson');

totalCorr = fastPearsonCorr(spikes1,spikes2);

%% shuffle spikes to calculate signal correlation

shuffledCorr=nan(nShuffle,1);

for ishuffle = 1:nShuffle
    shuffledSpikes = spikes1(randsample(nCounts,nCounts)); % shuffle spikes1
    
    % calculate correlation
    %     shuffledCorr(ishuffle) = ...
    %         (mean(shuffledSpikes.*spikes2) - mean(shuffledSpikes)*mean(spikes2))...
    %          /(std(shuffledSpikes).*std(spikes2));
    
    shuffledCorr(ishuffle) = fastPearsonCorr(shuffledSpikes,spikes2);
end


signalCorr = mean(shuffledCorr);

%% calculat noise corr

noiseCorr = totalCorr-signalCorr;

end
