% This file is meant to test the slicing of the reference file as well as the lining up of 
% records to this reference and their slicing afterwards
% Basically a test of the reference files
%
% christopher zaworski
% last edit : april 13, 2019
%
%


clc;close all;
disp('-----------reference_test.m------------')
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');

% path_folder = strcat('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/', folder)

references = {'/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/reference_files/031418_A0000B0000r27a.wav'};    


reference_file = references{1}
reference = audio_recordclass(reference_file)
reference.process_tracks();

time_offset = 0;
track_names = keys(reference.tracks) ;
track_data  = values(reference.tracks) ;
track_times = values(reference.track_times);


% reference_obj = save('031418_A0000B0000r27a.mat', reference_file, '-mat')
% reference.save_obj('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/reference_files/031418_A0000B0000r27a')


wave_name = reference.filename

% figure(1); hold on; grid on;
CSV_titles = {"wavefile","track_name", "RMS", "total_clicks", "common_clicks", "unique_clicks"};
CSV_MATRIX = cell(1,6);
CSV_MATRIX(1,:) = CSV_titles;
disp('before for loop')

size(CSV_MATRIX)
for i = (1:length(reference.tracks));
    % RMS_value = rms(track_data{i});
    RMS_value = -96;
%     plot(track_times{i}, track_data{i});
    % clicks = audio_clickdetect(reference.tracks('transition'), reference.fs);
    % [click_matrix, lagDiff] = audio_clickmatrix(clicks,clicks);

    % coarseness = 20; % number of samples we allow the mode to be off by when counting common clicks 
    % % [dt_row, dt_column] = find(click_matrix > (lagDiff - coarseness) & click_matrix < (lagDiff + coarseness));
    % clicks_total = size(click_matrix,1)
    % clicks_common = length(dt_row)
    % clicks_unique = size(clicks, 1) - length(dt_row)

    clicks_total = 20; 
    clicks_common = 10;
    clicks_unique = 3;
    size(CSV_MATRIX)
    size({wave_name, track_names{i}, RMS_value, clicks_total,clicks_common, clicks_unique})
    CSV_MATRIX = [CSV_MATRIX ; {wave_name, track_names{i}, RMS_value, clicks_total, clicks_common, clicks_unique}];

    % track_name = track_names{i}
    % size(CSV_MATRIX)
    % disp('CSV MATRIX')
    % CSV_MATRIX
    % CSV_MATRIX(end+1) = {wave_name, track_names, RMS_value, 'n/a', 'n/a', 'n/a'};
    % size(CSV_MATRIX)
end
disp('after for loop')
CSV_MATRIX
CSV_TABLE = cell2table(CSV_MATRIX)%,'VariableNames',CSV_titles)
writetable(CSV_TABLE,'A0000B0000_analysis.csv')