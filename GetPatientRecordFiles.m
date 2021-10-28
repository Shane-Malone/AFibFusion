% Script to obtain all records contained in the dataset folder
% Every patient record is contained in a folder named after the patient
% For example patient 608 has record files in "Datasets\p000608"

datasetFolder = 'C:\Users\shane\OneDrive - University College Dublin\Fifth Year\Project\Datasets';

patientFolders = dir(datasetFolder);
recordFile = fopen('patientRecords.txt', 'w');

for k = 3 : length(patientFolders) % First two folders are '.' and '..' so skip them
    patientName = patientFolders(k).name;
    fullFlderName = fullfile(patientFolders(k).folder, patientName);
    fprintf(1, 'Now reading patient: %s\n', patientName);
    
    patientFiles = dir(fullFlderName);
    for n = 3: length(patientFiles)
        baseFileName = patientFiles(n).name;
        % Obtain only the records named after the patient and ignore the
        % n.hea files
        if contains(baseFileName, patientName) && ~contains(baseFileName,"n.hea") && ~contains(baseFileName, ".atr")
            fprintf(1, '\tNow reading record: %s\n', baseFileName);
            fullFileName = "Datasets\"+patientName+"\"+baseFileName;
            fprintf(recordFile, '%s\n', fullFileName);
        end
    end
end
fclose(recordFile);