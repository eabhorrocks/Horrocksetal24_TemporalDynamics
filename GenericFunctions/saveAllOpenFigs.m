function saveAllOpenFigs(FolderName, FigNameprefix, FileType, renderer)

%FolderName = pwd; 
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = [FigNameprefix num2str(get(FigHandle, 'Number'))];
    set(0, 'CurrentFigure', FigHandle);
    switch FileType
        case '.fig'
            savefig(fullfile(FolderName, [FigName, '.fig']));
        otherwise
            print(FigHandle, [FigName], FileType, renderer)
    end
end

end