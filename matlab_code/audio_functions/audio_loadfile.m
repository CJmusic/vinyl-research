
%%% HOW TO HANDLE AUDIO ARRAYS in the form [data_ref, time_ref, fs_ref] 


%data_ref = load_audio[1];
%time_ref = load_audio[2];
%fs_ref = load_audio[3];


%data_ref[[LEFT,RIGHT]
%         [LEFT,RIGHT]
%            ...    
%         [LEFT,RIGHT]] 
%left_channel = data_ref[:,1];
%right_channel = data_ref[:,2];

%time_ref = [0.0 sec , 1/fs sec, 2/fs sec, ... , N];



DONT USE THIS 

%~~~~~~~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~%

%{ LOAD FILES %}
%This function takes a path and returns the data and time arrays, along with the sample rate

addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/');
function [data_ref, time_ref, fs_ref] = audio_loadfile(path_ref);
    [data_ref, fs_ref] = audioread(path_ref);
    time_ref = (0:length(data_ref)-1)/fs_ref;
end
%~ Loop through all the files and line them up with the reference
