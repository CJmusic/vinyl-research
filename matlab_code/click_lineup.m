%{
This file can also count the clicks in the record and see what clicks are common 
across all records (ones that are in the reference file as well), and ones that 
aren't. 

%}
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/260219_noisereferenceinst/');


addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/');
audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/'

AUDIO_FILES = {'one.wav','two.wav','three.wav', 'four.wav', 'five.wav'};
%AUDIO_FILES = {[fsinst_shorted, fsinst_shortedgained, gain10reference, gain10recordnoise, recordnoise, reference, system_noise, system_gained] };

path_ref = strcat(audio_dir,AUDIO_FILES{1});
path_ref
[data_ref, fs_ref] = audioread(path_ref);
time_ref = (0:length(data_ref)-1)/fs_ref;

% slice up the array into a much smaller more manageable chunck, 
% look for a local maximum 
% take 0.9 seconds worth of audio on either side of the max for one groove 

t_s = 5.8;
t_e = 6.2;

data_ref = data_ref(t_s*fs_ref:t_e*fs_ref,:);
%time_ref = time_ref(t_s*fs_ref:t_e*fs_ref);
time_ref = (0:length(data_ref)-1)/fs_ref; 

clf(figure(1));
%clf(figure(2));
%figure(1);
%grid on; hold on;
%plot(time_ref,data_ref,'g');


[ clicks_ref ] = audio_clickdetect(data_ref, fs_ref);
size(clicks_ref)
%figure(1);
%grid on; hold on;
%plot(time_ref,data_ref,'g');

%x = zeros(length(clicks_ref));
%figure(1);
%plot(clicks_ref/fs_ref,x, 'r.', 'MarkerSize', 20);

%figure(2);
%grid on; hold on;

for i = (1:length(AUDIO_FILES));         
    [data, time, fs] = audio_load(strcat(audio_dir,AUDIO_FILES{i}));
    data = data(t_s*fs:t_e*fs,:);
    time = (1:length(data))/fs;
    %[data, time] = audio_lineup(data, fs, time, data_ref);
    [ clicks ] = audio_clickdetect(data, fs);
    size(clicks)
    %lag = mode([[clicks] - [clicks_ref]]);
   
    %I want to find the closest matching value in the reference array, take the difference
    %then calculate the mode and correct for the lag, repeat this process a few times
    % I have two dimensions, amplitude and time. So I can look for the matched clicks that minimize
    % both the time difference and the amplitude difference

    cd_ref = diff(clicks_ref);
    cd_file = diff(clicks);
    values = [];
    amp_diffs = [];
    sam_diffs = [];
    for xi = (1:length(cd_ref));
        %[values, data_index] = min(abs(cd_ref(xi) - cd_file)); % subtract sample # of click arrays, take the absolute value, take the minimum
        [amp_diff, sam_diff] = min(abs(data_ref(cd_ref(xi)) - data(cd_file))); % subtract amplitude of click arrays, take abs, take min
        amp_diffs = [amp_diffs, amp_diff];
        sam_diffs = [sam_diffs, sam_diff]; 
    end
    
    amp_diffs
    sam_diffs



    figure(1);
    grid on; hold on;
    %plot(time,data, 'Color', [i/5,i/5,i/5] );
    x = zeros(length(clicks_ref));
    figure(1);
    plot(clicks_ref/fs,x, 'r.', 'MarkerSize', 20);

end

