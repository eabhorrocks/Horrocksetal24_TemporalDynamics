%% batch process sessions for analysis
dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';

sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

%%

saveSuffix = 'PSTH_fastSlow_v2.mat';
inputSuffix = 'basic.mat';

tic
parfor isession = 1:size(sessionTags,1)
    isession
    
    inputFileName = [sessionTags{isession,1},'_', sessionTags{isession,2},'_', inputSuffix];
    outputFileName = [sessionTags{isession,1},'_', sessionTags{isession,2},'_', saveSuffix];

      processSession_PSTH(inputFileName,outputFileName,dataDir) 
    % processSession_r2overTime(inputFileName,outputFileName,dataDir)
     %processSession_correlationAnalysis_v2(inputFileName,outputFileName,dataDir)
     %processSession_decoding(inputFileName,outputFileName,dataDir)
     %processSession_popAnalysis(inputFileName,outputFileName,dataDir)
    % processSession_FAoverTime(inputFileName,outputFileName,dataDir)

    %processSession_popAnalysis_fitSep(inputFileName,outputFileName,dataDir)




%     catch
%         disp(['failed on ', num2str(isession)]);
%         debug = 1;
%     end

end
toc