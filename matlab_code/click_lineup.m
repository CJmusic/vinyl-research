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

path_ref = strcat(audio_dir,AUDIO_FILES{2});
path_ref
[data_ref, fs_ref] = audioread(path_ref);
time_ref = (0:length(data_ref)-1)/fs_ref;

% slice up the array into a much smaller more manageable chunck, 
% look for a local maximum 
% take 0.9 seconds worth of audio on either side of the max for one groove 

t_s = 67.0;
t_e = 70.0;

data_ref = data_ref(t_s*fs_ref:t_e*fs_ref,:);
time_ref = (0:length(data_ref)-1)/fs_ref; 


[ clicks_ref ] = audio_clickdetect(data_ref, fs_ref);
size(clicks_ref)
%figure(1);
%grid on; hold on;
%plot(time_ref,data_ref,'g');
%for xi = 1:length(clicks_ref);
%    x1 = time_ref(clicks_ref(xi));
%    line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
%end
%
%x = zeros(length(clicks_ref));
%figure(1);
%plot(clicks_ref/fs_ref,x, 'r.', 'MarkerSize', 20);

%figure(2);
%grid on; hold on;

manual_clicks = []; %%sample numbers of the same click in each recording

clf(figure(1)); clf(figure(10));

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

    %this for loop attempts to look through the clicks in the reference and line them up with a click in the file
    diff_array = [];
    diag_diff = [];
    disp('for loop')
    for xi = (1:length(clicks_ref));%this makes an array with the distance between each click in the two
                                    % files
        diff_array = [diff_array; clicks - clicks_ref(xi)];
    end

    disp('mode')
    diff_array

    lagdiff = mode(diff_array(:))%size(diff_array))

    %% let's add some proper case handing for lagdiff > 0 and lagdiff < 0;
    if lagdiff > 0;
        % positive lagdiff means that the data array is delayed compared to the reference
        cdata = data(lagdiff + 1:end,:); 
        time = time - lagdiff/fs;%(1:length(data))*fs; 
    elseif lagdiff < 0; 
        % negative lagdiff means that the data array is ahead of the reference  
        cdata = data(1: end + lagdiff,:);   
        time = time - lagdiff/fs;%(1:length(data))*fs; 
    elseif lagdiff == 0; %then nothing needs to be corrected  
        cdata = data; %no lag 
        time = time - lagdiff/fs;%(1:length(data))*fs; 
    end

    %time = time - lagdiff/fs;%(1:length(data))*fs; 
    %cdata = data(lagdiff+1:end,:);
    size(cdata)
    size(data_ref)
    size_diff = length(data_ref) - length(cdata)
        
    if size_diff > 0; %% the following is to ensure the data arrays are the same length for calculating the coherence
                      %% probably this isn't necessary after the lag diff is properly handled 
        disp('>0'); 
        cdata_ref = data_ref(size_diff+1:end,:); 
    elseif size_diff < 0; 
        cdata = cdata(size_diff,:);
        cdata_ref = data_reff;
    else;
        cdata_ref = data_ref;
    end

    disp('sizes corrected')
    size(cdata)
    size(cdata_ref)

    figure(1);
    grid on; hold on;
    plot(time,data(:,1));%, 'Color', [i/5,i/5,i/5] );
    %plot(time_ref, data_ref(:,1), 'g');
    legend(AUDIO_FILES);
    %for xi = 1:length(clicks);
    %     x1 = time(clicks(xi));
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    %end
    title('Click detection: Amplitude vs. Time'); 
    xlabel('Time [s]');
    ylabel('Amplitude');

    figure(10); grid on; hold on;
    [amp_coh, freq_coh] = audio_mscohere(cdata_ref, cdata, fs);
    plot(freq_coh, amp_coh);
    legend(AUDIO_FILES);
    xlabel('frequency [Hz]')
    title('Coherence compared to first recording')
end


