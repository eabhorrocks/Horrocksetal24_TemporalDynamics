%% compare FA fit collectively or separataely

fileNames = {'FA_normal';...
    'FA_fitSep'};

plotTitles = {'fit together'; 'fit sepaarately'};

sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';

%% compare shared variance captured fittitng separately and together

for ifile = 1:numel(fileNames)
    ifile
    clear s
    s = struct;
    for isession = 1:5
        fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_', fileNames{ifile}, '.mat'];

        load(fullfile(dataDir,fname))

        s(isession).session = session;

        if ifile==1
            crit(ifile).SV(isession) = s(isession).session.s.SV;
            crit(ifile).qOpt(isession) = s(isession).session.s.qOpt;
            crit(ifile).q(isession) = s(isession).session.s.q;
        else
            crit(ifile).statSV(isession) = s(isession).session.stat.s.SV;
            crit(ifile).statqOpt(isession) = s(isession).session.stat.s.qOpt;
            crit(ifile).statq(isession) = s(isession).session.stat.s.q;

            crit(ifile).runSV(isession) = s(isession).session.run.s.SV;
            crit(ifile).runqOpt(isession) = s(isession).session.run.s.qOpt;
            crit(ifile).runq(isession) = s(isession).session.run.s.q;

            crit(ifile).SV(isession) = mean([crit(ifile).statSV(isession),crit(ifile).runSV(isession)]);
            crit(ifile).qOpt(isession) = mean([crit(ifile).statqOpt(isession),crit(ifile).runqOpt(isession)]);
            crit(ifile).q(isession) = mean([crit(ifile).statq(isession),crit(ifile).runq(isession)]);
        end
    end



end

compareSharedVariaance = cat(1,crit.SV)
[mean(compareSharedVariaance,2), sem(compareSharedVariaance,2)]


%% get velocity and speed for both types













