
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


if ismac() == true
    addpath('/Users/cz/Code/vinyl-research/matlab_code')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Wow')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
    addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    folder = ('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    RecordNumbers = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137_RecordNumbers.csv')

end 
if ispc() == true
    addpath('D:\Code\vinyl-research\matlab_code')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin')    
    addpath('E:\audio_files\A0137B0137')
    folder = ('E:\audio_files\A0137B0137\')
    RecordNumbers = readtable('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0137B0137\A0137B0137_RecordNumbers.csv')

end



pressingID = 'A0137B0137';


disp(['loading folder...:', folder])
% files = dir(strcat(folder,'*.wav'))
files = dir(fullfile(folder,'*.wav'))


% AudioTableHeaders = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};
AudioTableHeaders = {'pressing','record', 'side','track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicksa_L', 'commonclicksa_R','commonclicksb_L', 'commonclicksb_R', 'THD_L', 'THD_R', 'test_freq_L', 'wfreqspecamplitude_L', 'freqrms_L', 'WFrms_L', 'test_freq_R', 'wfreqspecamplitude_R', 'freqrms_R', 'WFrms_R',  'stereo_bleed'};


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
    files(i);
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
    
    file = strcat(files(i).folder,'/',filename);
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

    recordid = str2num(filename(1:end - 5));

    % pressid = RecordNumbers.pressing(strcmp(RecordNumbers.RecordID,recordid),:)
    pressing = RecordNumbers(RecordNumbers.RecordID == recordid, :);
    pressing = pressing.pressing{1};

    record = filename;
    side = filename(end-4);

    disp([strcat('...pressing:', pressing)])
    disp([strcat('...recordid:', recordid)])
    disp([strcat('...record:', record)])
    disp([strcat('...side:', side)])

    infoCell = {pressing, record, side};
    AudioOutput = RecordProcess(file);

    infoCell
    % cell2table([infoCell, AudioOutput], 'VariableNames', AudioTableHeaders)
    for j = (1:length(AudioOutput))
        AudioTable = [AudioTable; cell2table([infoCell, AudioOutput(j,:)], 'VariableNames', AudioTableHeaders)]
    end
    %append audio output to info cell array
    disp('SAVING CSV')
    writetable(AudioTable, strcat(folder,pressingID,'-AudioTable.csv'));
end


