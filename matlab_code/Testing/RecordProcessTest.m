close all; clear all; clc;


disp('~~~~~~~~~~~~TESTING RECORDPROCESS~~~~~~~~~~~~')

addpath('D:\Code\vinyl-research\matlab_code\')
% addpath('/Users/cz/Code/vinyl-research/matlab_code')


addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\')

addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
file = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\03141_A0000B0000r030b.wav'


% the reference file should be an input into record process (most likely a folder containing the wav file, and any relevent matrices as csv files) 



for i = (1:length(measurements))
    for j = (1:length(signal_names))
        % intTable = byN(strcmp(byN(:,9),signal_names{j}),:);
        int_stats = datastats(getData(byN, signal_names{j}, measurements{i}));

        intTable = cell2table({signal_names{j}, measurements{i},int_stats.max,int_stats.min,int_stats.mean,int_stats.median,int_stats.range,int_stats.std},'VariableNames', col_names);

        % intTable
        AudioStats = [AudioStats ; intTable];
    end
end
AudioStats

RecordProcess(file)

%% reference code from audio stats, appending to a table
for i = (1:length())
    for j = (1:length())
        % intTable = byN(strcmp(byN(:,9),signal_names{j}),:);
        int_stats = datastats(getData(byN, signal_names{j}, measurements{i}));

        intTable = cell2table({signal_names{j}, measurements{i},int_stats.max,int_stats.min,int_stats.mean,int_stats.median,int_stats.range,int_stats.std},'VariableNames', col_names);

        % intTable
        AudioStats = [AudioStats ; intTable];
    end
end
AudioStats