function indices = findContinousLogical(logicalVector, nElements)


spanLocs = bwlabel(logicalVector);   %identify contiguous logical true elements
spanLength = regionprops(spanLocs, 'area');  %length of each span of true elements
spanLength = [spanLength.Area];                         
goodSpans = find(spanLength>=nElements); % find any spans have >= than nElements
indices = find(ismember(spanLocs, goodSpans)); % get indices

end