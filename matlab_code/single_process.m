% This file is meant to process a single wav file for the 
% purpose of generating plots for presentations.
%
% christopher zaworski
% last edit : April 15, 2019
% 

clc;close all;
disp('-----------single_process.m------------')
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');

reference = audio_recordclass('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r28a.wav');
record = audio_recordclass('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r30a.wav');


addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')


% reference = audio_recordclass('/Volumes/AUDIOBANK/audio_files/A0000B0000/031418_A0000B0000r27a.wav');
% record = audio_recordclass('/Volumes/AUDIOBANK/audio_files/A0000B0000/03141_A0000B0000r27b.wav');


reference.process_tracks();
record.process_tracks();

if length(reference.clicks) > 0; 
    clicks_ref = reference.clicks;
else
    clicks_ref = audio_clickdetect(reference.tracks('transition'), reference.fs);
end


if length(reference.clicks) > 0; 
    clicks_tran = record.clicks;
else
    clicks_tran = audio_clickdetect(record.tracks('transition'), record.fs);
end
% clicks_tran = audio_clickdetect(record.tracks('transition'), record.fs);
data_refL = reference.tracks('transition');
data_refL = data_refL(:,1);

dataL = record.tracks('transition');
dataL = dataL(:,1);

figure(10); hold on; grid on; 
title('pre lineup')
plot(reference.track_times('transition'), data_refL); 
plot(record.track_times('transition'), dataL);


[click_matrix, lagdiff] = audio_clickmatrix(clicks_tran, clicks_ref);
record.lagdiff = lagdiff;
record.lagcorrect();
record.process_tracks();

data_refL = reference.tracks('transition');
data_refL = data_refL(:,1);

dataL = record.tracks('transition');
dataL = dataL(:,1);

figure(20); hold on; grid on; 
title('post lineup')
plot(reference.track_times('transition'), data_refL); 
plot(record.track_times('transition'), dataL);

data = record.tracks('transition');
data_ref = reference.tracks('transition');

coh_start = 10;
coh_end   = 70; 
rec_cohere = data(coh_start*record.fs:coh_end*record.fs,:);
ref_cohere = data_ref(coh_start*reference.fs:coh_end*reference.fs,:); 

figure(1); grid on;
n_sam = length(data_ref);
freq_fft = record.fs*(0:(n_sam/2))/n_sam;
data_fft = fft(data_ref)/n_sam;
data_fft = data_fft(1:n_sam/2+1);

plot(freq_fft, 20.0*log10(data_fft))  
set(gca, 'XScale', 'log');
title(strcat(reference.filename, "'s Spectrum"))
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  

figure(2); grid on; 
n_sam = length(data);
freq_fft = record.fs*(0:(n_sam/2))/n_sam;
data_fft = fft(data)/n_sam;
data_fft = data_fft(1:n_sam/2+1);

plot(freq_fft, 20.0*log10(data_fft))  
set(gca, 'XScale', 'log');
title(strcat(record.filename, "'s Spectrum"))
xlabel('Frequency (Hz)')
ylabel('Level (dB)') 


[ amp_coh, freq_coh ] = audio_mscohere(data, data_ref, reference.fs);

figure(3); grid on; 
plot(freq_coh,amp_coh(:,1))
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
% title('Coherences to Reference, Left Channel')
title(strcat('Coherences, Left Channel ', record.filename,' to ', reference.filename ))


figure(4); grid on; 
plot(freq_coh,amp_coh(:,2))
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
title(strcat('Coherences, Right Channel ', record.filename,' to ', reference.filename ))