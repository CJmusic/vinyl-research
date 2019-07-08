%addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/vinyl-research/matlab_code/audio_analysis.m');

addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/vinyl-research/matlab_code/');

filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.1.wav';

[data, time, fs] = load_audio(filename);


figure(1)
plot(time,data)


%%let's chop the audio around a click 
t_s = 6.94; %seconds
t_e = 6.95; %seconds 

s = find(t_s < time < t_e); 

data = data(s);
time = time(s);
%I need to detect the click, remove it, then median filter it
%
%classification of the click may need to be made based on its length
%
%
data_declicked = medfilt1(data,3);

figure(2)
plot(time,data_declicked)


function [data_ref, time_ref, fs_ref] = load_audio(path_ref)
    [data_ref, fs_ref] = audioread(path_ref);
    time_ref = (0:length(data_ref)-1)/fs_ref;
end
