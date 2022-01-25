annotationFileName = "C:\Users\shane\OneDrive - University College Dublin\Fifth Year\Project\AFibAnnotations\Annotations.csv";
annotations = readtable(annotationFileName);
annotationSize = size(annotations);

afibFilenames = annotations.('WF_FILENAME_1');

recordFile = fopen('patientRecords.txt', 'r');
Fs = 125;

% For each file create an N by 1 vector of annotations
% Create an N by 1 vector of sample locations for each annotaion

% Each dataFile represents a signal of varying length
dataFile = fgetl(recordFile);
while ischar(dataFile)
    
    fprintf(1, '\tNow writing annotations for record: %s\n', dataFile);
    
    for i=1:annotationSize(1)
        if contains(dataFile, afibFilenames(i))
            AFStatus = annotations{i, 'AFStatus'};
            Signalength = seconds(annotations{i, 'TotalTime_1'});
            
            annFilename = strrep(dataFile, '.hea', '');
            
            if AFStatus == 1
                OnsetTime = seconds(duration(annotations{i, 'OnsetTime_1'}));
                OffsetTime = seconds(duration(annotations{i, 'OffsetTime_1'}));
                startsample = OnsetTime*Fs;
                endsample = OffsetTime*Fs;
                
                annsamples = [startsample endsample];
                anntype = [];
                comments = {'A', 'N'};
                wrann(annFilename, 'atr', annsamples, anntype, 0, 0, 0, comments)
            else
                annsamples = 1;
                anntype = [];
                comments = {'N'};
                wrann(annFilename, 'atr', annsamples, anntype, 0, 0, 0, comments)
            end
            break
        end
    end
    dataFile = fgetl(recordFile);
end