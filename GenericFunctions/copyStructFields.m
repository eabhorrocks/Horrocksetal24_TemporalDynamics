function outStruct = copyStructFields(inStruct,outStruct)

if size(inStruct)~=size(outStruct)
    error('inStruct and outStruct have different sizes');
end

structFieldNames = fieldnames(inStruct);

for ifield = 1:numel(structFieldNames)
    newVals = num2cell(vertcat(inStruct.(structFieldNames{ifield})));
    [outStruct.(structFieldNames{ifield})] = newVals{:};
end
