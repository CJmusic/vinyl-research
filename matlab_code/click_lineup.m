
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
size(data_ref)
size(time_ref)

[val, idx] = max(data_ref(fs_ref*8.0:fs_ref*12.0,:)) %look for the max value in a small slice of the audio

data_ref_chop = data_ref(fs_ref*8.0 + idx(2)-0.9*fs_ref+1:fs_ref*8.0 + idx(2)+0.9*fs_ref,:);
time_ref_chop = time_ref(fs_ref*8.0 + idx(2)-0.9*fs_ref+1:fs_ref*8.0 + idx(2)+0.9*fs_ref);  


t_s = 5.8;
t_e = 6.2;

data_ref = data_ref(t_s*fs_ref:t_e*fs_ref,:);
%time_ref = time_ref(t_s*fs_ref:t_e*fs_ref);
time_ref = (0:length(data_ref)-1)/fs_ref; 

clf(figure(1));clf(figure(2));
figure(1);
grid on; hold on;
plot(time_ref,data_ref,'g');
audio_clickdetect(data_ref, fs_ref);

function clicks = audio_clickdetect(data, fs);
   time = (0:length(data)-1)/fs; 
   d_data = diff(data);%/diff(time);
   dd_data = diff(d_data);
   clicks = [];
   disp('for loop')
   for i = (1:length(data));
      if i + 1024 < length(data);
          threshold = rms(data(i:i+1024));
          %fft_wavelet = fft(wavelet);
          if d_data(i) > threshold;
                click = i
                i/fs  
                clicks = [clicks, click];   
          end
      end
   end 
   disp('for loop end')
   x = zeros(length(clicks));
   figure(1);
   plot(time, data); grid on; hold on;
   plot(clicks/fs,x, 'r.', 'MarkerSize', 20);
   figure(2); 
   subplot(2,1,1);
   plot(time(1:end-1),d_data); grid on; hold on;
   plot(clicks/fs,x,'r.', 'MarkerSize', 20);
   subplot(2,1,2); 
   plot(time(1:end-2),dd_data); grid on; hold on;
   plot(clicks/fs,x,'r.', 'MarkerSize', 20);
end 

%{
figure(1);
grid on; hold on;
plot(time_ref,data_ref,'g');

figure(2);
grid on; hold on;
plot(time_ref_chop,data_ref_chop,'g');

for i = (1:length(AUDIO_FILES));         
    [data, time, fs] = audio_load(strcat(audio_dir,AUDIO_FILES{i}));
    [data, time] = audio_lineup(data, fs, time, data_ref);
    %data = data(fs*8.0:fs*12.0,:); %slices the data array into a smaller section
    %time = time(fs*8.0:fs*12.0); 

    [val, idx] = max(data(fs*8.0:fs*12.0,:)); %look for the max value in a small slice of the audio
    
    data_chop = data(fs*8.0 + idx(2)-0.9*fs+1:fs*8.0 + idx(2)+0.9*fs,:); 
    time_chop = time(fs*8.0 + idx(2)-0.9*fs+1:fs*8.0 + idx(2)+0.9*fs);  

    figure(1);
    grid on; hold on;
    plot(time,data);
    
    figure(2);
    grid on; hold on;
    plot(time_chop,data_chop);

end
%}
