%This function takes a path and returns the data and time arrays, along with the sample rate

%%% HOW TO HANDLE AUDIO ARRAYS in the form [data, time, fs] 


%data = audio_load[1];
%time = audio_load[2];
%fs = audio_load[3];


%data[[LEFT,RIGHT]
%         [LEFT,RIGHT]
%            ...    
%         [LEFT,RIGHT]] 
%left_channel = data[:,1];
%right_channel = data[:,2];

%time = [0.0 sec , 1/fs sec, 2/fs sec, ... , N];


function [data, time, fs] = audio_load(path)
    addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/');
    [data, fs] = audioread(path);
    time = (0:length(data)-1)/fs;
    %data = data.';
end %audio_load
