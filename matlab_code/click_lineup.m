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

t_s = 5.0;
t_e = 10.0;

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

clf(figure(1));

for i = (1:2);%:length(AUDIO_FILES));         
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

    values = [];
    amp_diffs = [];
    sam_clicks = [];
    lag_diffs = [];
    %this for loop attempts to look through the clicks in the reference and line them up with a click in the file
    diff_array = [];
    diag_diff = [];
    disp('for loop')
    for xi = (1:length(clicks_ref));
        diff_array = [diff_array; clicks - clicks_ref(xi)];
        %[amp_diff, sam_click] = min(abs(data_ref(clicks_ref(xi)) - data(clicks))); % subtract amplitude of click arrays, take abs, take min
        %sam_click_ref = clicks_ref(xi);  % the sample that the click corresponds to in the ref data
        %lag_diff = sam_click - sam_click_ref;
        %lag_diffs = [lag_diffs, lag_diff];        
        %amp_diffs = [amp_diffs, amp_diff];
        %sam_clicks = [sam_clicks, sam_click]; 
    end
    size(diff_array)
    diff_array 
    for ii = (1:size(diff_array(1)-1));
        for j = (1:length(diff_array(ii,:)))
            if diff_array(ii,ii) - diff_array(ii+1, j) == 0;
                disp('ZERO')            
                if j-ii > 1;   
                 for k=(1:j-ii);
                  diff_array(ii+k,:) = [];  
              end
          end
      end 
        end
    end  
    diff_array
    %diag_diff
    %disp('TIME DELAY: ')
    %lag_diff/fs
    %length(diff_array)
    lagdiff = trace(diff_array)/length(diff_array);
    %lagdiff = clicks(1) - clicks_ref(i);
    %lagdiff = mode(diag_diff)
    %data = data(1:end-lagdiff);
    time = time + lagdiff/fs;%(1:length(data))*fs; 
    figure(1);
    grid on; hold on;
    plot(time,data(:,1));%, 'Color', [i/5,i/5,i/5] );

    for xi = 1:length(clicks);
         x1 = time(clicks(xi));
         line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    end
    title('Click detection: Amplitude vs. Time'); 
    xlabel('Time [s]');
    ylabel('Amplitude');
    %x = zeros(length(clicks_ref));
    %figure(1);
    %plot(clicks_ref/fs,x, 'r.', 'MarkerSize', 20);

end


