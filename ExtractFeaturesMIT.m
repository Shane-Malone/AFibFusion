%% Extract RR Intervals for each waveform
recordFile = fopen('files\RECORDS', 'r');
Fs = 250;
numSets = 1;
AF = cell(1); % Label for whole sequence
feature = cell(1);
windowSize = 20; % Beat intervals per window

lengthSegment = 496 ; % Length in seconds of signal to work on at a time
startTime = 0;
endTime = lengthSegment; % work on 2 minutes at a time
maxStartTime = 10000; % Max of 10000 seconds per record

valueset = ["AF","NonAF"];
% Each dataFile represents a signal of varying length
dataFile = fgetl(recordFile);
newRecord = true;
while ischar(dataFile)
    windowSample = 1;
    dataFile = strrep(dataFile, '.hea', '');
    if newRecord
        fprintf(1, 'Now extracting RR from record: %s\n', dataFile);
    end
    newRecord = false;

    startSample = startTime*Fs;
    endSample = endTime*Fs-1;

    try
        [signal,~,~] = rdsamp(strcat('files\', dataFile), [], endSample, startSample, 1); % Read ECG signal
        [ann,~,~,~,~,comments] = rdann(strcat('files\', dataFile), 'atr', [], [], [], []);  % Read entire annotations
    catch
        warning('Cannot read from record. Moving to next record.')
        dataFile = fgetl(recordFile);
        startTime = 0;
        endTime = lengthSegment;
        newRecord = true;
        continue
    end

    [RRIntervalSet, secLocs] = RRFinder(signal, Fs); % Read entire set of intervals and corresponding samples
    
    if isempty(RRIntervalSet) % If no beats found move to next record
        dataFile = fgetl(recordFile);
        continue
    end
    
    maxNumWindows = 1000; % Extract max of 10000 windows of size windowSize of intervals
    if sum(contains(comments, '(AFIB')) > 0
        AFibInComments = true;
    end
    
    while windowSample < length(RRIntervalSet)-windowSize && windowSample < maxNumWindows*windowSize

        % Create labels for each RR Interval
        AF(numSets) = {0};
        AFFound = false;
        if AFibInComments
            for sampleOfBeat = secLocs(windowSample:windowSample+windowSize-1) % Where each beat is in the segment being analysed
                globalSample = sampleOfBeat + startSample; % Where each beat is in relation to whole signal
                % Find valid annotation
                for currAnn = 1:length(ann)
                    if globalSample > ann(currAnn)
                        if currAnn == length(ann)
                            break
                        elseif globalSample < ann(currAnn + 1)
                            break
                        end
                    end
                end
                if strcmp(comments{currAnn}, '(AFIB')
                    AF(numSets) = {1};
                    break;
                end
            end
        end

        feature(numSets) = {RRIntervalSet(windowSample:windowSample+windowSize-1)}; % Pass this window of samples
        numSets = numSets+1;
        windowSample = windowSample + windowSize;
    end
    
    startTime = startTime + lengthSegment;
    endTime = startTime + lengthSegment;
%     if startTime >= maxStartTime
%         % Get next data file
%         dataFile = fgetl(recordFile);
%         startTime = 0;
%         endTime = lengthSegment;
%         newRecord = true;
%     end
    
end

save('FeatureSetMIT20Beats', 'feature', 'AF')