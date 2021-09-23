close all; clear all;



if ismac() == true
    addpath('/Users/cz/Code/vinyl-research/matlab_code')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/A0137B0137')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
    addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
end 
if ispc() == true
    addpath('D:\Code\vinyl-research\matlab_code')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin')    
    addpath('E:\audio_files\A0137B0137')
end

addpath('/Users/cz/Code/vinyl-research/matlab_code/Wow')


% file = 'd:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/lacquer_recordings/Dec 20 - Test A Part one .wav'
% file = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/lacquer_recordings/lacquerpartone.wav'
% offset = 12.6;

% lacquerprocess(file, offset)

folder = '/Volumes/AUDIOBANK/audio_files/lacquer_recordings';

pressingID = 'lacquer';


disp(['loading folder...:', folder])
% files = dir(strcat(folder,'*.wav'))
files = dir(fullfile(folder,'*.wav'))


% AudioTableHeaders = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};
AudioTableHeaders = {'pressing','record', 'side','track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicksa_L', 'commonclicksa_R','commonclicksb_L', 'commonclicksb_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};


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
    filename = files(i).name;
    files(i);
    disp(['opening file...:', filename])

    
    if ismember(filename, AudioTable.record)
        disp('record already processed...')
        continue
    end
    
    file = strcat(files(i).folder,'/',filename);
    date_recorded = 0;
    pressing = 0;
    top_stamper = 0;
    top_hits = 0;
    bottom_stamper = 0;
    bottom_hits = 0;
    record = 0;
    side  = 0;
    track  = 0;
    
    recordid = 0;
    record = filename;
    pressing = 'lacquer'

    disp([strcat('...pressing:', pressing)])
    disp([strcat('...recordid:', recordid)])
    disp([strcat('...record:', record)])
    disp([strcat('...side:', side)])

    infoCell = {pressing, record, side};
    if i == 1;
        AudioOutput = LacquerProcess(file, 19.3);
    end
    if i == 2; 
        AudioOutput = LacquerProcess(file, 0.1);
    end

    AudioOutput
    infoCell
    length(AudioOutput)
    sz = size(AudioOutput)
    loop = sz(1)
    % cell2table([infoCell, AudioOutput], 'VariableNames', AudioTableHeaders)
    for j = (1:loop)
        AudioTable = [AudioTable; cell2table([infoCell, AudioOutput(j,:)], 'VariableNames', AudioTableHeaders)];
    end
    %append audio output to info cell array
    disp('SAVING CSV')
    writetable(AudioTable, strcat(folder, '/', pressingID,'-AudioTable.csv'));
end