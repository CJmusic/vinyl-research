% This file looks at the specific grooves on a single
% recordings and does the measurements we need for that
% 
% christopher zaworski
%
% last edited : april 25, 2019

clc; clear all; close all;
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
disp('-----------groove_process.m------------')
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');

addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')


% reference = audio_recordclass('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r28a.wav');
% record = audio_recordclass('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r30a.wav');
tstart = 10;
tend = 25;

[data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/030419_r26fivetrials/one.wav');
% [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/042619_lacquerApostsweepv.wav');

data = data(tstart*fs:tend*fs,:);


time = linspace(0,(length(data)-1)/fs,length(data));
rotation_speed = 33.33333;%45;
T = 60/rotation_speed; %this is the length of one groove segment
num_segs = (floor(length(data)/fs/T))
n_sam = round(T*fs)
time_seg = time(1:n_sam);
seg_array = []; %need to 


for ng = 1:num_segs
    seg_array(:,:,ng) = data(1+(ng-1)*n_sam:ng*n_sam,:);
end

figure(1);
grid on; hold on;

for ng = 1:num_segs-1
    %%% coherences in the left channel 
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