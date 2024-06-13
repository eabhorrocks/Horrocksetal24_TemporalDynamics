function units = genUnitStructFromAllenData(unitInformation,unitIndex,spikeTimes)

% inputs are:
% unitInfo: (1x1) struct with (1xnUnit) fields describing different details of units
% unitIndex: nUnit vector, index of unitInfo fields
% spikeTimes: (1x1) struct with nUnit fieldNames of unit IDs and vectors of
% spiketimes

unitIDs = fieldnames(spikeTimes);
units = struct();
unitInfoNames = fieldnames(unitInformation);
nUnits = numel(unitIDs);

for iunit = 1:nUnits
    units(iunit).ID = str2num(erase(unitIDs{iunit}, 'c_'));
    units(iunit).spiketimes = spikeTimes.(unitIDs{iunit});
    
    % spikeTimes and unitInfo don't have same order
    idx = find(unitIndex==units(iunit).ID); 
    for ifield = 1:numel(unitInfoNames)
        units(iunit).(unitInfoNames{ifield}) = unitInformation.(unitInfoNames{ifield})(idx);
    end
end



end