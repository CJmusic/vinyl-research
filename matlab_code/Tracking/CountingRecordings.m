

%~~~~~~~~~~~  COUNTING A0000B0000 BEGINS ~~~~~~~~~~~~%
 
addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')
folder = ('/Volumes/AUDIOBANK/audio_files/A0000B0000/')


timestamps = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000timestamps.csv');

head(timestamps)

recording = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/Spreadsheets/Tracking_Sheets/RecordingTracking.xlsx');

head(recording)



for i = (1:height(recording))
    % recording{i,3}
    i
    recording.ASide{i}
    recording.ASide{i} = 'n';
    recording.BSide{i} = 'n';
    recording.AUDIOBANK{i} = 'n';
    recording.OneDrive{i} = 'n';
    
end

head(recording)


disp(['loading folder...:', folder])
% files = dir(strcat(folder,'*.wav'))
files = dir(fullfile(folder,'*.wav'))

% for i = (1:height(recording))
    
% end


for i = (1:length(files)) %%loop through records
    filename = files(i).name;
    % disp(filename)
    recordid = str2num(filename(19:21));
    disp(recordid)
    side = filename(22);
    disp(side)


    % index = ismember(recording.RECORDID, str2num(filename))

    % index = str2num(recordid)
    if strcmp(side, 'a');
        recording.ASide{recordid} = 'y';
    end
    if strcmp(side, 'b');
        recording.BSide{recordid} = 'y';
    end


end

head(recording)
writetable(recording, 'A0000B0000recordingtracking.csv')

%~~~~~~~~~~~  COUNTING A0000B0000 ENDS ~~~~~~~~~~~~%

%~~~~~~~~~~~  COUNTING A0137B0137 BEGINS ~~~~~~~~~~~~%
 
addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
folder = ('/Volumes/AUDIOBANK/audio_files/A0137B0137/')


% timestamps = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000timestamps.csv');

% head(timestamps)

recording = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/Spreadsheets/Tracking_Sheets/RecordingTracking.xlsx','Sheet', 2);

head(recording)



for i = (1:height(recording))
    % recording{i,3}
    i
    recording.ASide{i}
    recording.ASide{i} = 'n';
    recording.BSide{i} = 'n';
    recording.AUDIOBANK{i} = 'n';
    % recording.OneDrive{i} = 'n';
    
end

head(recording)


disp(['loading folder...:', folder])
% files = dir(strcat(folder,'*.wav'))
files = dir(fullfile(folder,'*.wav'))

% for i = (1:height(recording))
    
% end


for i = (1:length(files)) %%loop through records
    filename = files(i).name;
    disp(filename)
    recordid = (filename(1:3));
    disp(recordid)
    side = filename(4);
    disp(side)


    % index = ismember(recording.RECORDID, str2num(filename))
    % index = find((str2num(recording.RecordID, recordid)));
    % index

    % index = str2num(recordid)
    recording.ASide(recording.RecordID == str2num(recordid))

    if strcmp(side, 'a');
        recording.ASide(recording.RecordID == str2num(recordid)) = {'y'};
    end
    if strcmp(side, 'b');
        recording.BSide(recording.RecordID == str2num(recordid)) = {'y'};
    end


end

writetable(recording, 'A0137B0137recordingtracking.csv')

%~~~~~~~~~~~  COUNTING A0137B0137 ENDS ~~~~~~~~~~~~%