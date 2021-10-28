% This file looks at the specific grooves on a single
% recordings and does the measurements we need for that
% 
% christopher zaworski
%
% last edited : april 25, 2019

clc; close all; clear all; 
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
disp('-----------groove_process.m------------')
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');

addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/Common/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/from_John')


% [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/030419_r26fivetrials/one.wav');
% [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/042619_lacquerApostsweepv.wav');
% [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/lacquer+pioneer/052319_pioneer.wav');

% tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/-03b1558.270.wav')

% tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0000B0000/081219-A0000B0000r093b1555.849.wav')

tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0000B0000/071919-A0000B0000r038b1559.173.wav')

data = tracks('quiet');
fs = 96000;

% 1 khz
tstart = 6;
tend = 21;
data = data(tstart*fs:tend*fs,:);
time = linspace(0,(length(data)-1)/fs,length(data));
rotation_speed = 33.33333;%45;
% n_sam = 172800;
n_sam = 1.8*96000;
num_segs = (floor(length(data)/n_sam))
time_seg = time(1:n_sam);
seg_array = []; 

for ng = 1:num_segs
    seg_array(:,:,ng) = data(1+(ng-1)*n_sam:ng*n_sam,:);
end

for ng = 1:num_segs-1
    %%% coherences in the left channel 
    disp('NEW LOOP~~~~~~~~~~~~~~~~~~~~~~~~~')
    %%********* PRINT THE LAGDIFF BETWEEN GROOVES ************
    if ng > 1;
        groovediff = audio_corrlineup(seg_array(:,1,ng), seg_array(:,1,ng-1));
        seg_array(:,2,ng-1) = circshift(seg_array(:,2,ng-1), groovediff);
    end

    %%%%%
    % coh_LR = audio_mscohere(seg_array(:,1,ng), seg_array(:,2,ng),fs); 

    [coh_nextL, freq_coh] = audio_mscohere(seg_array(:,1,ng), seg_array(:,1,ng+1), fs); % previous left channel with next left channel
    [coh_firstL, ~] = audio_mscohere(seg_array(:,1,1), seg_array(:,1,ng+1), fs); % first left channel with next left channel

    %%% coherences in the right channel    
    [coh_nextR, ~] = audio_mscohere(seg_array(:,2,ng), seg_array(:,2,ng+1), fs);% previous right channel with next right channel
    [coh_firstR, ~] = audio_mscohere(seg_array(:,2,1), seg_array(:,2,ng+1), fs);% first right channel with next right channel

    [coh_nextLR, ~] = audio_mscohere(seg_array(:,1,ng), seg_array(:,2,ng+1), fs); % previous left channel with next right channel
    [coh_nextRL, ~] = audio_mscohere(seg_array(:,2,ng), seg_array(:,1,ng+1), fs); % previous right channel with next left channel


    [coh_LR, ~] = audio_mscohere(seg_array(:,1,ng), seg_array(:,2,ng), fs);% same groove, left and right channel

    coh_firstL = pwroctsmooth_singlesided(coh_firstL, 0.33);
    coh_firstR = pwroctsmooth_singlesided(coh_firstR, 0.33);
    coh_nextL = pwroctsmooth_singlesided(coh_nextL, 0.33);
    coh_nextR = pwroctsmooth_singlesided(coh_nextR, 0.33);
    coh_nextLR = pwroctsmooth_singlesided(coh_nextLR, 0.33);
    coh_nextRL = pwroctsmooth_singlesided(coh_nextRL, 0.33);
    coh_LR = pwroctsmooth_singlesided(coh_LR, 0.33);


    figure(1); hold on;
    subplot(2,1,1)
    hold on; grid on;
    plot(time_seg,seg_array(:,1,ng), 'DisplayName',['segment',num2str(ng)])
    title('left channel') 
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')

    subplot(2,1,2)
    hold on; grid on;
    plot(time_seg,seg_array(:,2,ng), 'DisplayName',['segment',num2str(ng)]) 
    title('right channel') 

    % plot(time_seg,seg_array(:,1,ng), 'DisplayName',['segment',num2str(ng)], 'Color', [1 - 1.0*ng/100 ,0,1.0*ng/100,0.5]);\
    title('Groove segments')
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('time (s)')
    ylabel('Amplitude')

    figure(2); hold on;
    subplot(2,1,1)
    hold on;
    plot(freq_coh, coh_nextL, 'DisplayName',['segment',num2str(ng)])
    title('left channel')
    xlabel('Frequency (Hz)')
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    subplot(2,1,2)
    hold on;
    plot(freq_coh, coh_nextR, 'DisplayName',['segment',num2str(ng)])
    title('right channel')
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    % title(strcat('Coherences next groove Left Channel'))


    figure(3); hold on;
    plot(freq_coh, coh_nextLR, 'DisplayName',['segment',num2str(ng)])
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    title(strcat('Coherences Left and adjacent Right Channel'))


    figure(4); hold on;
    plot(freq_coh, coh_nextRL, 'DisplayName',['segment',num2str(ng)])
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    title(strcat('Coherences Right and next Left Channel'))



    figure(5); hold on;
    subplot(2,1,1)
    hold on;
    plot(freq_coh, coh_firstL, 'DisplayName',['segment',num2str(ng)])
    title('left channel')
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    subplot(2,1,2)
    hold on;
    plot(freq_coh, coh_firstR, 'DisplayName',['segment',num2str(ng)])
    title('right channel')
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')


    % title(strcat('Coherences first groove Left Channel'))

    figure(6); hold on;
    plot(freq_coh, coh_LR, 'DisplayName',['segment',num2str(ng)])
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    title(strcat('Coherence in noise between Left and Right Channel'))


    figure(7); hold on; %% plot the spectrum of each groove
    subplot(2,1,1)
    hold on;
    n_sam = length(seg_array(:,1,ng));
    freq_fft = fs*(0:(n_sam/2))/n_sam;
    data_fft = fft(seg_array(:,1,ng))/n_sam;
    data_fft = data_fft(1:n_sam/2+1);

    data_fft = pwroctsmooth_singlesided(data_fft,0.33);

    plot(freq_fft, 20.0*log10(data_fft)) 
    grid on; 
    set(gca, 'XScale', 'log');
    title('left channel')    
    xlabel('Frequency (Hz)')
    ylabel('Level (dB)')  
    subplot(2,1,2)
    hold on;
    n_sam = length(seg_array(:,2,ng));
    freq_fft = fs*(0:(n_sam/2))/n_sam;
    data_fft = fft(seg_array(:,2,ng))/n_sam;
    data_fft = data_fft(1:n_sam/2+1);

    data_fft = pwroctsmooth_singlesided(data_fft,0.33);

    plot(freq_fft, 20.0*log10(data_fft)) 
    grid on; 
    set(gca, 'XScale', 'log');
    title('right channel')    
    xlabel('Frequency (Hz)')
    ylabel('Level (dB)')  



    [cor_first, lags_first] = xcorr(seg_array(:,1,1), seg_array(:,1,ng));
    [cor_next, lags_next] = xcorr(seg_array(:,1,ng), seg_array(:,1,ng+1));
    
    figure(8); hold on;
    plot(lags_first, cor_first);
    grid on; 
    title('xcorr with first groove')
    xlabel('lag')
    ylabel('correlation')

    figure(9); hold on;
    plot(lags_next, cor_next);
    grid on; 
    title('xcorr with next groove')
    xlabel('lag')
    ylabel('correlation')
end

saveas(figure(1),'groovesegments.png')
saveas(figure(2),'coherencenext.png')
saveas(figure(3),'coherenceLR.png')
saveas(figure(4),'coherenceRL.png')
saveas(figure(5), 'coherencefirst.png')
saveas(figure(6), 'coherencestereo.png')
saveas(figure(7),'groovespectrum.png')