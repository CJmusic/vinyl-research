% This file calculated the coherence between multiple files 
%
% christopher zaworski
% last edit : march 31, 2019
%
%
clear all; clc;close all;
disp('------------freqcohere.m---------------')

addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/260219_noisereferenceinst/');

%addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/');
%audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/'
%
%AUDIO_FILES = {'one.wav','two.wav','three.wav', 'four.wav', 'five.wav'};

addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/');
audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/';

AUDIO_FILES = {'Bcorrelation_test_1.wav','Bcorrelation_test_2.wav','Bcorrelation_test_3.wav'};

addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');
audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/';
AUDIO_FILES = {'031418_A0000B0000r29a-cutsilence.wav','031418_A0000B0000r27a-cutsilence.wav','031418_A0000B0000r28a-cutsilence.wav'}
path_ref = strcat(audio_dir,AUDIO_FILES{1});










%%% The code below is legacy, and does not use the implementation of the record class
%{
path_ref
[data_ref, fs_ref] = audioread(path_ref);
%data_ref = data_ref(:,1);
%data_ref = data_ref(1.0*fs_ref:6.0*fs_ref,:);
time_ref = (0:length(data_ref)-1)/fs_ref;

ref_coh = data_ref;%(length(data_ref)/4:length(data_ref)/4+2^20,:);

clicks_ref = audio_clickdetect(data_ref, fs_ref);
cdata_ref  = data_ref;
%cdata_ref = data_ref(7.0*fs_ref : 15*fs_ref);


% the loop below is based on manual click detection
manual_clicks = [7.203497085, 11.50687415, 16.336197884] % these are the timestamps for John's Bcorrelation files
%cdata_manual_ref = data_ref(manual_clicks(1)*fs_ref:manual_clicks(1)*fs_ref + 15.0*fs_ref);
%ctime_ref = (0:length(cdata_manual_ref))/fs_ref;
for i = (1:length(AUDIO_FILES));
    strcat(audio_dir,AUDIO_FILES{i})
    [data, time, fs] = audio_load(strcat(audio_dir,AUDIO_FILES{i}));
    %data = data(:,1);
    data = data(1.0*fs:6.0*fs,:);

    clicks = audio_clickdetect(data, fs);
    [click_matrix, lagDiff] = audio_clickmatrix(clicks, clicks_ref);
    lagDiff

    % if lagDiff > 0;
    %     cdata = data(lagDiff:end);
    % end
    % if lagDiff < 0;
    %     cdata = data(1:lagDiff);
    % end
    % if lagDiff == 0;
    %     cdata = data;
    % end
    cdata = circshift(data, lagDiff);
    length(cdata)
    ctime = (1:length(cdata))/fs;
    length(ctime)
    figure(1); hold on; grid on;
    plot(ctime, cdata)
    title('Automatic lining up of Audio via clickmatrix')
    legend(AUDIO_FILES)

    cdata = cdata(1*fs : 5*fs,:);
    cdata_ref = data_ref(1*fs_ref : 5*fs_ref,:);

    [ amp_coh, freq_coh ] = audio_mscohere(cdata_ref, cdata, fs_ref);
    figure(2); grid on; hold on;
    plot(freq_coh,amp_coh(:,1))
    set(gca, 'XScale', 'log');
    %axis([0,fs_ref/2,0,1])
    xlabel('frequency [Hz]')
    title('Groove Coherences, Left')
    legend(AUDIO_FILES)

    figure(3); grid on; hold on;
    plot(freq_coh,amp_coh(:,2))
    set(gca, 'XScale', 'log');
    %axis([0,fs_ref/2,0,1])
    xlabel('frequency [Hz]')
    title('Groove Coherences, Right')
    legend(AUDIO_FILES)

    figure(4);
    time = (1:length(data))*fs;
    plot(time, data);
    for xi = 1:length(clicks);
         x1 = time(clicks(xi));
         line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
     end

    %% code below is to manually line up the clicks
    % cdata_manual = data(manual_clicks(i)*fs:manual_clicks(i)*fs + 15.0*fs);
    % ctime_manual = (1:length(cdata_manual))/fs;
    % [ amp_coh, freq_coh ] = audio_mscohere(cdata_manual_ref, cdata_manual, fs_ref);
    % nfft=2^14;
    % freq_coh = ([0:nfft/2])*fs_ref/nfft;
    % amp_coh = mscohere(cdata_manual_ref,cdata_manual,hanning(nfft),[],nfft,fs_ref);

    % figure(5); hold on; grid on;
    % plot(ctime_manual,cdata_manual)
    % title('Manually Lining up Audio')
    % legend(AUDIO_FILES)

    % figure(4); grid on; hold on;
    % plot(freq_coh,amp_coh(:,1))
    % set(gca, 'XScale', 'log');
    % %axis([0,fs_ref/2,0,1])
    % xlabel('frequency [Hz]')
    % title('Groove Coherences, Manual')
    % legend(AUDIO_FILES)

end
%}