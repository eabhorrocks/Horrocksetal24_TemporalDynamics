function centreOfMass = calcCOM(cellArray,xvals, method)


if nargin == 3
    if strcmp(method, 'mean')
        xSums = cellfun(@mean, cellArray);
    end
else
    xSums = cellfun(@sum, cellArray); % default
end

    
    M = sum(xSums(:));

    centreOfMass = sum((xvals.*xSums))/M;