function trials = getRunAndPupilInfo(trials, runInfo, pupilInfo, options)

if ~exist('options','var'),                options=struct;                end
% timing options
if ~isfield(options,'preTime'),            options.preTime=0.2;            end
if ~isfield(options,'postTime'),           options.postTime=0.8;           end
if ~isfield(options,'sampleSpacing'),      options.sampleSpacing=0.01;           end

% smoothing options
if ~isfield(options,'smoothMethod'),         options.smoothMethod='gaussian';  end
if ~isfield(options,'smoothWindow'),       options.smoothWindow=175;           end


% resample run and pupil traces to 10ms bin size and smooth run speed
sampleSpacing = 0.01; %10ms
queryTimes = runInfo.runMidTimes(1):sampleSpacing:runInfo.runMidTimes(end);
runInfo.interpTime = queryTimes';
runInfo.rawSpeedInterp = interp1(runInfo.runMidTimes, runInfo.runSpeed, queryTimes)';
windowSize_bin = options.smoothWindow/(sampleSpacing*1000);
runInfo.smthSpeedInterp = smoothdata(runInfo.rawSpeedInterp,options.smoothMethod,windowSize_bin);

if isstruct(pupilInfo)
sampleSpacing = 0.01; %10ms
queryTimes = pupilInfo.pupilTime(1):sampleSpacing:pupilInfo.pupilTime(end);
pupilInfo.interpTime = queryTimes';
pupilInfo.interpPupilDilation = interp1(pupilInfo.pupilTime, pupilInfo.pupilDilation, queryTimes)';
pupilInfo.interpPupil_x = interp1(pupilInfo.pupilTime, pupilInfo.pupil_x, queryTimes)';
pupilInfo.interpPupil_y = interp1(pupilInfo.pupilTime, pupilInfo.pupil_y, queryTimes)';

pupilInfo.interpPupilDilation = interp1(pupilInfo.pupilTime, pupilInfo.pupilDilation, queryTimes)';
pupilInfo.interpPupilDilation = interp1(pupilInfo.pupilTime, pupilInfo.pupilDilation, queryTimes)';

end

for itrial = 1:numel(trials)

    startrunidx = find(runInfo.interpTime<trials(itrial).start_time-options.preTime,1,'last');
    stoprunidx = find(runInfo.interpTime>trials(itrial).stop_time+options.postTime,1,'first');

    trials(itrial).runTrace = runInfo.smthSpeedInterp(startrunidx:stoprunidx);
    trials(itrial).runTime = runInfo.interpTime(startrunidx:stoprunidx);
    trials(itrial).meanRunSpeed = mean(trials(itrial).runTrace);
    
    if isstruct(pupilInfo) % some sessions don't have eye tracking
      
        startpupilidx = find(pupilInfo.interpTime<trials(itrial).start_time-options.preTime,1,'last');
       stoppupilidx = find(pupilInfo.interpTime>trials(itrial).stop_time+options.postTime,1,'first');

        trials(itrial).pupilDilationTrace = pupilInfo.interpPupilDilation(startpupilidx:stoppupilidx);
        trials(itrial).meanPupilDilation = nanmean(trials(itrial).pupilDilationTrace);
        
        trials(itrial).pupilTrace = [pupilInfo.interpPupil_x(startpupilidx:stoppupilidx),...
                                        pupilInfo.interpPupil_y(startpupilidx:stoppupilidx)];
                                    
%         trials(itrial).pupilArea = pupilInfo.pupilArea(startpupilidx:stoppupilidx);
%         trials(itrial).meanPupilArea = nanmean(trials(itrial).pupilArea);
        trials(itrial).pupilTime = pupilInfo.interpTime(startpupilidx:stoppupilidx);

    else 
        trials(itrial).pupilDilationTrace = NaN;
        trials(itrial).meanPupilDilation = NaN;
        trials(itrial).pupilTrace = NaN;
%         trials(itrial).meanPupilArea = NaN;
        trials(itrial).pupilTime = NaN;
    end
end

end