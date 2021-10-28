%% Extract RR Intervals for each waveform
recordFile = fopen('patientRecords.txt', 'r');
pulse_prop_delay = 64; %pulse propagation delay is -75 for this signal
Fs = 125;
numSets = 1;
AFperSample = cell(1); % per sample label
AF = cell(1); % Label for whole sequence
feature = cell(1);

startTime = 60;
endTime = startTime + 240;
startSample = Fs*startTime;
endSample = Fs*endTime - 1;

% Each dataFile represents a signal of varying length
dataFile = fgetl(recordFile);
while ischar(dataFile)
    dataFile = strrep(dataFile, '.hea', '');
    fprintf(1, '\tNow extracting RR from record: %s\n', dataFile);
    
    [signal,~,~] = rdsamp(dataFile, [], endSample, startSample, 3); % Read entire ECG signal
    signal(:,1) = signal(:,1)-mean(signal(:,1));
    [ann,~,~,~,~,comments] = rdann(dataFile, 'atr', [], [], [], []);  % Read entire annotations
    
    [RRIntervalSet, secLocs, ~, ~, ~, ~] = RRFinder(signal, Fs, pulse_prop_delay);
    
    % Create labels for each RR Interval
    AF_sequence = zeros(1, length(RRIntervalSet));
    AF(numSets) = {0};
    for i=1:length(RRIntervalSet)
        sampleOfBeat = secLocs(i) + startSample;
        [~, annIndex] = ismember(sampleOfBeat, ann);
        if annIndex > 0
            if comments(annIndex) == "(AFIB"
                AF_sequence(i) = 1;
                AF(numSets) = {1};
            end
        end
    end
    if isempty(RRIntervalSet)
        dataFile = fgetl(recordFile);
        continue
    end
    AFperSample(numSets) = {AF_sequence};

    feature(numSets) = {RRIntervalSet};
    numSets = numSets+1;
    
    % Get next data file
    dataFile = fgetl(recordFile);
end

save('FeatureSet', 'feature', 'AF', 'AFperSample')
