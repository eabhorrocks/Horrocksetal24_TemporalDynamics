

nGood_stat = histcounts(stat.unit(statGoodVals),'BinEdges',0.5:1:10000);
nGood_run = histcounts(stat.unit(runGoodVals),'BinEdges',0.5:1:10000);

idx = find(nGood_stat>=3 & nGood_run>=3);

% nGood = histcounts(stat.unit(bothGoodVals),'BinEdges',0.5:1:10000);
% idx = find(nGood>=4);


%%

idx = [291, 303, 367, 384, 424, 494, 558, 736, 744, 936, 975, 1078];
% tthese are ignoring sesh 3!

for iunit = 1:numel(idx)
    figure
    for ispeed = 1:6
        subplot(2,3,ispeed), hold on

          shadedErrorBar(1:200,mean(stat.psthVec(stat.unit==idx(iunit) & stat.speedVec==ispeed,:),1),...
         sem(stat.psthVec(stat.unit==idx(iunit) & stat.speedVec==ispeed,:),1),...
         'lineProps', 'k')
     shadedErrorBar(1:200,mean(run.psthVec(run.unit==idx(iunit) & run.speedVec==ispeed,:),1),...
         sem(run.psthVec(run.unit==idx(iunit) & run.speedVec==ispeed,:),1),...
         'lineProps', 'r')

     mygca(ispeed) = gca;
     title(ispeed)
    end

    title(idx(iunit))

    yl = cell2mat(get(mygca, 'Ylim'));
    ylnew = [min(yl(:,1)) max(yl(:,2))];
    set(mygca, 'Ylim', ylnew)
%      pause
%      close
end