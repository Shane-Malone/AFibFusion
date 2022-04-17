%% Extract RR Intervals for each waveform
% Will add white gaussian noise to evvery fifth window of 496 seconds
SignalLocs = readtable('ECG_PPG_SignalLocations.csv');
Fs = 125;
numSets = 1;
AF = cell(1); % Label for whole sequence
feature = cell(1);
windowSize = 20; % Beat intervals per window

lengthSegment = 496 ; % Length in seconds of signal to work on at a time
startTime = 0;
endTime = lengthSegment; % work on 2 minutes at a time
maxStartTime = 120000; % Max of 10000 seconds per record

valueset = ["AF","NonAF"];
% Each dataFile represents a signal of varying length
recordNum = 1;
newRecord = true;

windowNum = 1;

%for recordNum = [9, 14]
while recordNum < size(SignalLocs, 1)
    if recordNum ~= 9 && recordNum ~= 14
        recordNum = recordNum + 1;
        continue
    end
    windowSample = 1;
    dataFile = SignalLocs{recordNum, 'Record'};
    dataFile = dataFile{1};
    ecgLoc = SignalLocs{recordNum, 'ECG'};
    ppgLoc = SignalLocs{recordNum, 'PPG'};
    pulseDelay = SignalLocs{recordNum, 'Delay'};
    if newRecord
        fprintf(1, 'Now extracting RR from record: %s\n', dataFile);
        windowNum = 1;
    else
        windowNum = windowNum + 1;
        if windowNum >= 5
            windowNum = 1;
        end
    end
    newRecord = false;

    startSample = startTime*Fs;
    endSample = endTime*Fs-1;

    try
        [signal,~,~] = rdsamp(dataFile, [], endSample, startSample, 1); % Read ECG signal
        [ann,~,~,~,~,comments] = rdann(dataFile, 'atr', [], [], [], []);  % Read entire annotations
    catch
        warning('Cannot read from record. Moving to next record.')
        recordNum = recordNum + 1;
        startTime = 0;
        endTime = lengthSegment;
        newRecord = true;
        continue
    end

    ecgSig = signal(:,ecgLoc);
    ppgSig = signal(:,ppgLoc);

    SNR = -15;
    if windowNum == 1
        ecgSig = awgn(ecgSig, SNR);
        %ppgSig = awgn(ppgSig, SNR);
    end

    [RRIntervalSet, secLocs] = ECG_PPG_RRFinder(ecgSig, ppgSig, Fs, pulseDelay); % Read entire set of intervals and corresponding samples
    
    if isempty(RRIntervalSet) % If no beats found move to next record
        startTime = startTime + lengthSegment;
        endTime = startTime + lengthSegment;
        continue
    end
    
    maxNumWindows = 1000; % Extract max of 10000 windows of size windowSize of intervals
    if sum(contains(comments, 'A')) > 0
        AFibInComments = true;
    else
        AFibInComments = false;
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
                if strcmp(comments{currAnn}, 'A')
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
    if startTime >= maxStartTime
        % Get next data file
        recordNum = recordNum + 1;
        startTime = 0;
        endTime = lengthSegment;
        newRecord = true;
    end
    
end

save('FeatureSets2/FusedFeatureSetMIMIC20Beats_15dB', 'feature', 'AF')