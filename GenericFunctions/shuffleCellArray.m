function shuf = shuffleCellArray(cellArray,nPerms)

allCounts = vertcat(cellArray{:});
shuf(nPerms).array = [];

for iperm = 1:nPerms
    shufCounts = allCounts(randperm(numel(allCounts)));
    idx = 1;
    for icond = 1:numel(cellArray)
        nEle = numel(cellArray{icond});
        shuf(iperm).array{icond} = shufCounts(idx:((idx+nEle)-1));
        idx = idx+ nEle;
    end
    shuf(iperm).array = reshape(shuf(iperm).array,size(cellArray));
end

