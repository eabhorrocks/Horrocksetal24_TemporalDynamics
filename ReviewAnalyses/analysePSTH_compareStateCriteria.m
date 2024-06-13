%% compare PSTH responses using diff state criteria

fileNames = {'PSTHdata_wCellTypes.mat';...
             'PSTHdata_strict.mat';...
             'PSTHdata_cp.mat'};

plotTitles = {'original'; 'strict'; 'changepoints'};

%% loop through filenaames
for ifile = 1:numel(fileNames)
    clear stat run
% load data
load(fileNames{ifile})

% get reliable responses

rangeThresh = 3;
reliThrsh =-1.645;
zThresh = 3.29;

statGoodVals =stat.reliVec<=reliThrsh & stat.rangeFRVec>=rangeThresh & stat.peakAbsZVec>=zThresh;
runGoodVals = run.reliVec<=reliThrsh & run.rangeFRVec>=rangeThresh & run.peakAbsZVec>=zThresh;
bothGoodVals = statGoodVals & runGoodVals;

% plot average responses
figure(22), 
subplot(1,3,ifile),
hold on
shadedErrorBar(1:200, mean(stat.psthVec(statGoodVals,:),1), sem(stat.psthVec(statGoodVals,:),1))
shadedErrorBar(1:200, mean(run.psthVec(runGoodVals,:),1), sem(run.psthVec(runGoodVals,:),1), 'lineProps', {'Color', 'r'})
title(plotTitles{ifile})
ax = gca; ax.XTick=0:20:200; ax.XTickLabel = -0.2:0.2:1.8;
ylabel('Firing rate (Hz)')
xlabel('time')
defaultAxesProperties(gca, false)


% plot sustaianedness index
metric_name = 'susIndexVec';
%figure, hold on
statVals = stat.(metric_name);
runVals = run.(metric_name);
diffVals{ifile} =  runVals(bothGoodVals)-statVals(bothGoodVals);
%boxplot([statVals(bothGoodVals), runVals(bothGoodVals)])
%title(plotTitles{ifile})

end

%% plot delta sustainedness index for each

figure, hold on
distributionPlot(diffVals,'histOpt',1,'globalNorm',3,'color',[.7 .7 .7],'showMM',6,'divFactor',3)
plot([0 4],[0 0],'k:')
ylabel('delta Sustainedness index (run-stat)')
xlabel('State criteria')
ax = gca; ax.XTick = 1:3; ax.XTickLabel = plotTitles;
