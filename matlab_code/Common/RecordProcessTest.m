
close all; clear all; clc;


disp('~~~~~~~~~~~~TESTING RECORDPROCESS~~~~~~~~~~~~')

% % ~~~~ WINDOWS ~~~~ %
if ispc() == true 
    addpath('D:\Code\vinyl-research\matlab_code\')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
    addpath('D:\Code\vinyl-research\matlab_code\Wow\')
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\')
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000')
    addpath('E:\audio_files\A0000B0000\')


    % folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0137B0137\';
    % folder = 'E:\audio_files\A0000B0000\';
    folder = ('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\macvsacer\')
end
% % ~~~~ WINDOWS ENDS ~~~~ %



% ~~~~ MAC ~~~~ %
if ismac() == true
    addpath('/Users/cz/Code/vinyl-research/matlab_code')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')
    % file = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/003141_A0000B0000r30a.wav'

    % folder = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/';
    % addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    % folder = ('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    addpath('/Volumes/AUDIOBANK/audio_files/duplicaterecordingtest/')
    folder = ('/Volumes/AUDIOBANK/audio_files/duplicaterecordingtest/')
end

% ~~~~ MAC ENDS ~~~~ %

pressingID = 'A0137B0137'


disp(['loading folder...:', folder])
% files = dir(strcat(folder,'*.wav'))
files = dir(fullfile(folder,'*.wav'))


% AudioTableHeaders = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};
AudioTableHeaders = {'record', 'side','track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicksa_L', 'commonclicksa_R','commonclicksb_L', 'commonclicksb_R',  'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};


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

    filename = files(i).name;
    files(i)

    disp(['opening file...:', filename])

    
    if ismember(filename, AudioTable.record)
        disp('record already processed...')
        continue
    end
    
    file = strcat(files(i).folder,'/',filename)
    % file = strcat(files(i).folder,'\',filename)

    date_recorded = 0;
    pressing = 0;
    top_stamper = 0;
    top_hits = 0;
    bottom_stamper = 0;
    bottom_hits = 0;
    record = 0;
    side  = 0;
    track  = 0;
    

    record = filename;
    side = 'a';
    infoCell = {record, side};
    AudioOutput = RecordProcess(file);

    
    infoCell
    AudioOutput
    AudioTableHeaders


    for j = (1:length(AudioOutput))
        AudioTable = [AudioTable; cell2table([infoCell, AudioOutput(j,:)], 'VariableNames', AudioTableHeaders)];
    end
    %append audio output to info cell array
    disp('SAVING CSV')
    writetable(AudioTable, strcat(folder,pressingID,'-AudioTable.csv'));
end


