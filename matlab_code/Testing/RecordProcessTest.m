close all; clear all; clc;


disp('~~~~~~~~~~~~TESTING RECORDPROCESS~~~~~~~~~~~~')

% ~~~~ WINDOWS ~~~~ %
% addpath('D:\Code\vinyl-research\matlab_code\')
% addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
% addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\')
% file = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\03141_A0000B0000r030b.wav'
% folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\';

% ~~~~ MAC ~~~~ %
addpath('/Users/cz/Code/vinyl-research/matlab_code')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')
file = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/003141_A0000B0000r30a.wav'

folder = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/';




disp(['loading folder...:', folder])
files = dir(strcat(folder,'*.wav'))

AudioTableHeaders = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};


% check if there is already a csv file to append to 
try
    disp('trying...')
    strcat(folder,'A0000B0000-AudioTable.csv')
    AudioTable = readtable(strcat(folder,pressingID,'-AudioTable.csv'));
catch
    disp('csv file not found, creating one...')
    AudioTable  = cell2table(cell(0,length(AudioTableHeaders)), 'VariableNames', AudioTableHeaders);
end

for i = (1:length(files)) %%loop through records
    AudioTable
    AudioTable.record
    i
    filename = files(i).name

    disp(['opening file...:', filename])
    filename(19:21)
    AudioTable.record
    if ismember(filename(19:21), AudioTable.record)
        disp('record already processed...')
        continue
    end
    file = strcat(files(i).folder,'/',files(i).name);

    date_recorded = 0;
    pressing = 0;
    top_stamper = 0;
    top_hits = 0;
    bottom_stamper = 0;
    bottom_hits = 0;
    record = 0;
    side  = 0;
    track  = 0;
    
    %STRIP RELEVANT INFO FROM NAME 
    '031418_A0000B0000r27a.wav'
     123456789012345678901

    date_recorded = (filename(1:6))
    % date_recorded = date_recorded{1};
    record = filename(19:20)
    record = str2num(filename(19:20))

    top_stamper = filename(8)
    pressing = filename(8:17)
    top_hits = str2num(filename(9:12)) + record;
    bottom_stamper = filename(13);
    bottom_hits = str2num(filename(14:17)) + record;
    side = filename(21);

    top_hits = num2str(top_hits);
    bottom_hits = num2str(bottom_hits);
    record = num2str(record);

    disp([strcat('...date_recorded:', date_recorded)])
    disp([strcat('...pressing:', pressing)])
    disp([strcat('...record:', record)])
    disp([strcat('...top_stamper:', top_stamper)])
    disp([strcat('...top_hits:', top_hits)])
    disp([strcat('...bottom_stamper:', bottom_stamper)])
    disp([strcat('...bottom_hits:', bottom_hits)])
    disp([strcat('...side:', side)])

    infoCell = {date_recorded, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, record, side}; 
    disp('PROCESSING OUTPUT')
    AudioOutput = RecordProcess(file)
    size(AudioOutput)
    for j = (1:length(AudioOutput))
        disp('looping through audio output')
        j
        AudioOutput(j,1)
        AudioTable = [AudioTable; cell2table([infoCell, AudioOutput(j,:)], 'VariableNames', AudioTableHeaders)]
    end
    %append audio output to info cell array
    % infoCell
    % AudioOutput
    % AudioTableHeaders
    % AudioTable = [AudioTable; cell2table([infoCell, AudioOutput], 'VariableNames', AudioTableHeaders)]
    writetable(AudioTable, strcat(folder,pressingID,'-AudioTable.csv'));

    
end



%% reference code from audio stats, appending to a table
% for i = (1:length(Records))
%         % RecordProcess = ;
%         intTable = cell2table({},'VariableNames', col_names);

%         % intTable
%         output = RecordProcess(file);
%         output
%         AudioStats = [AudioStats ; intTable];
% end
% AudioStats

% Audio = cell2table(csvdata(2:end,:),'VariableNames',csvdata(1,:))
