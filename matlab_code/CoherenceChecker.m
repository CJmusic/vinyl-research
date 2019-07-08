%% This file will measure the coherences in grooves of a record 



clc; %close all; clear all; 
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
disp('-----------groove_process.m------------')
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');

addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')

% [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/030419_r26fivetrials/one.wav');
% [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/042619_lacquerApostsweepv.wav');
[data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/lacquer+pioneer/052319_pioneer.wav');



%signal_names = {'leadin','1kHz', '10kHz', '100Hz', 'freqsweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition', '1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2', '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout'};
%timestamps = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
%lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz

% 1 khz
tstart = 3;
tend = 10;

% quiet track on record
%tstart = 180;
%tend = 246;
data = data(tstart*fs:tend*fs,:);

time = linspace(0,(length(data)-1)/fs,length(data));
rotation_speed = 33.33333;%45;
% T = 60/rotation_speed; %this is the length of one groove segment
n_sam = 172800;
% num_segs = (floor(length(data)/fs/T))
num_segs = (floor(length(data)/n_sam))
% n_sam = round(T*fs)
time_seg = time(1:n_sam);
seg_array = []; %need to 


for ng = 1:num_segs
    seg_array(:,:,ng) = data(1+(ng-1)*n_sam:ng*n_sam,:);
end
for ng = 1:num_segs-1

    %%% coherences in the left channel 
    disp('NEW LOOP~~~~~~~~~~~~~~~~~~~~~~~~~')
    %%********* PRINT THE LAGDIFF BETWEEN GROOVES ************
    if ng > 1;
        groovediff = audio_corrlineup(seg_array(:,1,ng), seg_array(:,1,ng-1), fs)
        seg_array(:,2,ng-1)= circshift(seg_array(:,2,ng-1), groovediff);
    end

    %%%%%
    [coh_nextL, freq_coh] = audio_mscohere(seg_array(:,1,ng), seg_array(:,1,ng+1), fs);
    [coh_firstL, ~] = audio_mscohere(seg_array(:,1,1), seg_array(:,1,ng+1), fs);

    %%% coherences in the right channel    
    [coh_nextR, ~] = audio_mscohere(seg_array(:,2,ng), seg_array(:,2,ng+1), fs);
    [coh_firstR, ~] = audio_mscohere(seg_array(:,2,1), seg_array(:,2,ng+1), fs);


    figure(1); hold on;
    plot(time_seg,seg_array(:,1,ng), 'DisplayName',['segment',num2str(ng)]) 
    % plot(time_seg,seg_array(:,1,ng), 'DisplayName',['segment',num2str(ng)], 'Color', [1 - 1.0*ng/100 ,0,1.0*ng/100,0.5]);
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('time (s)')
    ylabel('Amplitude')

    figure(2); hold on;
    plot(freq_coh, coh_nextL, 'DisplayName',['segment',num2str(ng)])
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    title(strcat('Coherences next groove Left Channel'))


    figure(3); hold on;
    plot(freq_coh, coh_nextR, 'DisplayName',['segment',num2str(ng)])
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    title(strcat('Coherences next groove Right Channel'))


    figure(4); hold on;
    plot(freq_coh, coh_firstL, 'DisplayName',['segment',num2str(ng)])
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    title(strcat('Coherences first groove Left Channel'))

    figure(5); hold on;
    plot(freq_coh, coh_firstR, 'DisplayName',['segment',num2str(ng)])
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    title(strcat('Coherences first groove Right Channel'))


    figure(6); hold on; %% plot the spectrum of each groove
    n_sam = length(seg_array(:,1,ng));
    freq_fft = fs*(0:(n_sam/2))/n_sam;
    data_fft = fft(seg_array(:,1,ng))/n_sam;
    data_fft = data_fft(1:n_sam/2+1);

    plot(freq_fft, 20.0*log10(data_fft)) 
    grid on; 
    set(gca, 'XScale', 'log');
    title(strcat("Groove's Spectrum"))
    xlabel('Frequency (Hz)')
    ylabel('Level (dB)')  

    [cor_first, lags_first] = xcorr(seg_array(:,1,1), seg_array(:,1,ng));
    [cor_next, lags_next] = xcorr(seg_array(:,1,ng), seg_array(:,1,ng+1));
    
    figure(7); hold on;
    plot(lags_first, cor_first);
    grid on; 
    title('xcorr with first groove')
    xlabel('lag')
    ylabel('correlation')

    figure(8); hold on;
    plot(lags_next, cor_next);
    grid on; 
    title('xcorr with next groove')
    xlabel('lag')
    ylabel('correlation')
end

figure(1)
legend()

figure(2)
legend()

figure(3)
legend()

figure(4)
legend()

figure(5)
legend()

figure(6)
legend()

figure(7)
legend()

figure(8)
legend()