function outTable = flexibleTableRead(filename)

%% Notes
% this function is designed primarily for csv files produced by bonsai
% datatype is assumed to be numeric except for the last column which is 
% assumed to be a timestmap

% filename = 'X:\ibn-vision\DATA\SUBJECTS\M23002\SDTraining\230313\TrialParams2023-03-13T16_32_24.csv';

fid = fopen(filename);
linenum=2;
C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
C = strsplit(C{1}{:},',');
nVars = numel(C);

intialTable = readtable(filename);
nVars = width(intialTable);
nTrials = height(intialTable);

varNames=[];
varType=[];

for ivar = 1:nVars-1 %last column is datetime
    exampleCell = C{ivar};
    varNames{ivar} = exampleCell(isstrprop(exampleCell,'alpha'));
    varType = cat(2,varType, "double");
end

varNames{nVars} = 'Time';
varType = cat(2,varType, 'datetime');

% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", nVars);

% Specify range and delimiter
opts.DataLines = [2, inf];
opts.Delimiter = ",";
opts.VariableTypes = varType;
opts.VariableNames = string(varNames);
%opts = setvaropts(opts, ["StartSession", "VarName2", "VarName3", "VarName4", "VarName5"], "TrimNonNumeric", true);
opts = setvaropts(opts, string(varNames(1:end-1)), "TrimNonNumeric", true);

outTable = readtable(filename, opts);