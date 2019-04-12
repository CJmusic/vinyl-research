% This file is meant to test the slicing of the reference file as well as the lining up of records to this reference and their slicing afterwards
%
% christopher zaworski
% last edit : april 11, 2019
%
%


disp('----------reference_test.m-------------')
clc;close all;
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');

% path_folder = strcat('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/', folder)

references = {'/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav'};    


reference_file = references{1}
reference = audio_recordclass(reference_file)
reference.process_tracks();

time_offset = 0;
track_names = keys(reference.tracks) ;
track_data  = values(reference.tracks) ;
track_times = values(reference.track_times);

figure(1); hold on; grid on;
for i = (1:length(reference.tracks));
    % track = track_data{i}.';%'
    % time = (time_offset:time_offset + length(track)-1)/reference.fs;%
    % if i > 1;
    %     time_offset = time_offset + time(end);
    % end
    % time(1)
    % time(end)
    % size(time)
    % time_offset = time_offset + length(time);
    % time = time + time_offset;
    % time_offset
    % figure(i)
    plot(track_times{i}, track_data{i});
    % time = [];
end