
close all; clear all; clc;


disp('~~~~~~~~~~~~TESTING RECORDPROCESS~~~~~~~~~~~~')

% % ~~~~ WINDOWS ~~~~ %
% addpath('D:\Code\vinyl-research\matlab_code\')
% addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
% addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\')
% addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000')
% addpath('E:\audio_files\A0000B0000\')


% folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0137B0137\';
% folder = 'E:\audio_files\A0000B0000\';
% folder = ('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\macvsacer\')

% % ~~~~ WINDOWS ENDS ~~~~ %


<<<<<<< HEAD

% ~~~~ MAC ~~~~ %
addpath('/Users/cz/Code/vinyl-research/matlab_code')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')
addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/Vinyl_Project/data/A0137B0137')
% file = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/003141_A0000B0000r30a.wav'

% folder = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/';
addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
folder = ('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
% addpath('/Volumes/AUDIOBANK/audio_files/duplicaterecordingtest/')
% folder = ('/Volumes/AUDIOBANK/audio_files/duplicaterecordingtest/')
=======
if ismac() == true
    addpath('/Users/cz/Code/vinyl-research/matlab_code')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
    addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    folder = ('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
end 
if ispc() == true
    addpath('D:\Code\vinyl-research\matlab_code')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin')    
    addpath('E:\audio_files\A0137B0137')
    folder = ('E:\audio_files\A0137B0137')
end
>>>>>>> c8a099e08159e5c3d8aa38fd2530f894932356e6



pressingID = 'A0137B0137'


disp(['loading folder...:', folder])
% files = dir(strcat(folder,'*.wav'))
files = dir(fullfile(folder,'*.wav'))


% AudioTableHeaders = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};
<<<<<<< HEAD
AudioTableHeaders = {'pressing','record', 'side','track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};
RecordNumbers = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137_RecordNumbers.csv')
=======
AudioTableHeaders = {'pressing','record', 'side','track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicksa_L', 'commonclicksa_R','commonclicksb_L', 'commonclicksb_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};


>>>>>>> c8a099e08159e5c3d8aa38fd2530f894932356e6
% check if there is already a csv file to append to 
try
    disp('trying...')
    strcat(folder,strcat(pressingID,'-AudioTable.csv'))
    AudioTable = readtable(strcat(folder,strcat(pressingID,'-AudioTable.csv')))
catch
    disp('csv file not found, creating one...')
    AudioTable  = cell2table(cell(0,length(AudioTableHeaders)), 'VariableNames', AudioTableHeaders);
end
for i = (1:length(files)) %%loop through records
    % AudioTable
    % AudioTable.record
    % i
    filename = files(i).name;
    files(i)
    % filename = filename(3:end);
    disp(['opening file...:', filename])
    % filename(19:21)
    % AudioTable.record
    % filename(19:21)
    % num2str(AudioTable.record)
    % ismember(filename(19:21), num2str(AudioTable.record))
    % disp('ismember')
    % ismember(str2num(filename(19:21)), (AudioTable.record))
    % ismember(str2num(filename(22)), (AudioTable.side))
    
    if ismember(filename, AudioTable.record)
        disp('record already processed...')
        continue
    end
    
    file = strcat(files(i).folder,'/',filename)
    % file = strcat(files(i).folder,'\',filename)

    % pressid = 0; 
    date_recorded = 0;
    pressing = 0;
    top_stamper = 0;
    top_hits = 0;
    bottom_stamper = 0;
    bottom_hits = 0;
    record = 0;
    side  = 0;
    track  = 0;
    
    % PressingNumber 
    % RecordID 
    % pressing 
    % RecordNumber 

    recordid = str2num(filename(1:3));

    % pressid = RecordNumbers.pressing(strcmp(RecordNumbers.RecordID,recordid),:)
    pressing = RecordNumbers(RecordNumbers.RecordID == recordid, :)
    pressing = pressing.pressing{1}

    record = filename;
    % pressid = filename(1:end-6); % verify these 2 lines MAR 13 2020 
    side = filename(end-4);
    % pressing = filename(1:strfind(filename,string(sscanf(string(filename),'%*[^0123456789]%d')))-1)



    % top_stamper = filename(8);
    % pressing = filename(8:17);
    % top_hits = str2num(filename(9:12)) + record;
    % bottom_stamper = filename(13);
    % bottom_hits = str2num(filename(14:17)) + record;
    % side = filename(22);

    % top_hits = num2str(top_hits);
    % bottom_hits = num2str(bottom_hits);
    % record = num2str(record);

    % disp([strcat('...date_recorded:', date_recorded)])
    % disp([strcat('...pressing:', pressing)])
    % disp([strcat('...record:', record)])
    % disp([strcat('...top_stamper:', top_stamper)])
    % disp([strcat('...top_hits:', top_hits)])
    % disp([strcat('...bottom_stamper:', bottom_stamper)])
    % disp([strcat('...bottom_hits:', bottom_hits)])
    % disp([strcat('...side:', side)])

    % infoCell = {str2num(date_recorded), pressing, top_stamper, str2num(top_hits), bottom_stamper, str2num(bottom_hits), str2num(record), side}; 
    infoCell = {pressing,record, side};
    AudioOutput = RecordProcess(file);
    disp('AUDIOTABLE')
    AudioTable
    % cell2table([infoCell, AudioOutput], 'VariableNames', AudioTableHeaders)
    for j = (1:length(AudioOutput))
        AudioTable = [AudioTable; cell2table([infoCell, AudioOutput(j,:)], 'VariableNames', AudioTableHeaders)];
    end
    %append audio output to info cell array
    disp('SAVING CSV')
    writetable(AudioTable, strcat(folder,pressingID,'-AudioTable.csv'));
end


