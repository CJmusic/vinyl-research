close all; clear all; clc;


disp('~~~~~~~~~~~~TESTING RECORDPROCESS~~~~~~~~~~~~')

addpath('D:\Code\vinyl-research\matlab_code\')
% addpath('/Users/cz/Code/vinyl-research/matlab_code')


addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\')

addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
file = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\03141_A0000B0000r030b.wav'

folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\';
disp(['loading folder...:', folder])
files = dir(strcat(folder,'*.wav'));

AudioTableHeaders = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};

for i = (1:length(files))
    filename = files(i).name;
    disp(['opening file...:', filename])

    % T.record
    if isempty(AudioTable.record)
    elseif ismember(str2num(filename(19:21)), AudioTable.record)
        disp('record already processed...')
        continue
    end
    % disp('T{7}')
    % T
    % % if isempty(T{8,:})
    % if ismember(str2num(filename(19:21)), [T{:,7}])
    %     disp('record already processed...')
    %     continue
    % end

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
    date_recorded = (filename(1:6));
    % date_recorded = date_recorded{1};
    record = str2num(filename(19:21));

    top_stamper = filename(8);
    pressing = filename(8:16);
    top_hits = str2num(filename(9:12)) + record;
    bottom_stamper = filename(13);
    bottom_hits = str2num(filename(14:17)) + record;
    side = filename(22);

    top_hits = num2str(top_hits);
    bottom_hits = num2str(bottom_hits);
    record = num2str(record);

    disp([strcat('...date_recorded:', date_recorded)])
    disp([strcat('...pressing:', pressing)])
    disp([strcat('...:record', record)])
    disp([strcat('...top_stamper:', top_stamper)])
    disp([strcat('...top_hits:', top_hits)])
    disp([strcat('...bottom_stamper:', bottom_stamper)])
    disp([strcat('...bottom_hits:', bottom_hits)])
    disp([strcat('...side:', side)])

    infoCell = {date_recorded, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, record, side}; 
    
    AudioOutput = RecordProcess(file);
    %append audio output to info cell array
    AudioTable = cell2table({infoCell, AudioOutput}, 'VariableNames', AudioTableHeaders)
    
end



%% reference code from audio stats, appending to a table
for i = (1:length(Records))
        % RecordProcess = ;
        intTable = cell2table({},'VariableNames', col_names);

        % intTable
        output = RecordProcess(file);
        output
        AudioStats = [AudioStats ; intTable];
end
AudioStats

T = cell2table(csvdata(2:end,:),'VariableNames',csvdata(1,:))
writetable(AudioTable,'AudioTable.csv');
