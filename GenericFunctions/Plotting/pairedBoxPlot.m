function [bp, a] = pairedBoxPlot(nestedCellArray, subGroupSpacing, colorCellArray, defaultOpts)

%nestedCellArray = {{statta.nRespVec},{runta.nRespVec}};
%colorCellArray = {{cols}, {cols*0.7}};
%subGroupSpacing = [-0.15, +0.15];

nGroups = numel(nestedCellArray);
nSubGroups = cellfun(@numel, nestedCellArray);

if numel(unique(nSubGroups))~=1
    error('nSubGroups must be equal for each group')
end

nSubGroups = nSubGroups(1);

% build vector for box plot
allVals = [];

for isubgroup = 1:nSubGroups
    for igroup = 1:nGroups
        subgroup(isubgroup).group(igroup).vals = nestedCellArray{igroup}{isubgroup};
        subgroup(isubgroup).group(igroup).idx = ...
            (isubgroup+subGroupSpacing(igroup))*ones(size(nestedCellArray{igroup}{isubgroup}));
        subgroup(isubgroup).group(igroup).vals = subgroup(isubgroup).group(igroup).vals(:);
        subgroup(isubgroup).group(igroup).idx = subgroup(isubgroup).group(igroup).idx(:);
        
    end
    subgroup(isubgroup).allVals = vertcat(subgroup(isubgroup).group(:).vals);
    subgroup(isubgroup).allIdx = vertcat(subgroup(isubgroup).group(:).idx);
end

%% generate the plot

if nargin>3 
    if defaultOpts

bp = boxplot(vertcat(subgroup.allVals),vertcat(subgroup.allIdx),...
    'positions',vertcat(subgroup.allIdx),...
    'BoxStyle', 'outline', 'widths', 0.1,'outliersize', 4);
ax = gca;
a = ax.Children.Children;
for isubgroup = 1:nSubGroups
    for ielement = 1:numel(a)
        for igroup = 1:nGroups
            if a(ielement).XData(1) == isubgroup+subGroupSpacing(igroup)
                a(ielement).Color = colorCellArray{igroup}{1}(isubgroup,:);
                a(ielement).MarkerEdgeColor = colorCellArray{igroup}{1}(isubgroup,:);
            end
        end
    end
end
    end
end

if nargin < 4
    bp = boxplot(vertcat(subgroup.allVals),vertcat(subgroup.allIdx),...
    'positions',vertcat(subgroup.allIdx),...
    'BoxStyle', 'filled', 'widths', 0.1,'outliersize', 4, 'whisker', 1000);
ax = gca;
a = ax.Children.Children; %get(get(gca,'children'),'children');



for isubgroup = 1:nSubGroups
    for ielement = 1:numel(a)
        for igroup = 1:nGroups
            if a(ielement).XData(1) == isubgroup+subGroupSpacing(igroup)
                a(ielement).Color = colorCellArray{igroup}{1}(isubgroup,:);
                a(ielement).MarkerEdgeColor = colorCellArray{igroup}{1}(isubgroup,:);
            end
        end
    end
end

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'w');
set(lines, 'LineWidth', 2)
outs = findobj(gcf, 'type', 'line', 'Tag', 'Outliers');
set(outs, 'Marker', 'o');
t = get(a,'tag');   % List the names of all the objects
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',8); % Set width
ax.TickDir = 'out'; box off
%ax.XTick = 1:nSubGroups; 

end
