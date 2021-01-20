clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')

folder = '/Volumes/AUDIOBANK/audio_files/samerecordtest/';

files = dir(fullfile(folder,'*.wav'));

[data1, fs] = audioread('/Volumes/AUDIOBANK/audio_files/samerecordtest/1.wav');
[data2, fs] = audioread('/Volumes/AUDIOBANK/audio_files/samerecordtest/2.wav');
[data3, fs] = audioread('/Volumes/AUDIOBANK/audio_files/samerecordtest/3.wav');
[data4, fs] = audioread('/Volumes/AUDIOBANK/audio_files/samerecordtest/4.wav');
[data5, fs] = audioread('/Volumes/AUDIOBANK/audio_files/samerecordtest/5.wav');

% data1 = data1(8*fs:end,:);
% data2 = data2(8*fs:end,:);
% data3 = data3(8*fs:end,:);
% data4 = data4(8*fs:end,:);
% data5 = data5(8*fs:end,:);
figure(100)

[acor_L,lags_L] = xcorr(data1(:,1),data2(:,1));
[M_L,I_L] = max(abs(acor_L));
lagdiff12 = lags_L(I_L);
plot(lags_L, acor_L); hold on;

[acor_L,lags_L] = xcorr(data1(:,1),data3(:,1));
[M_L,I_L] = max(abs(acor_L));
lagdiff13 = lags_L(I_L);
plot(lags_L, acor_L); hold on;

[acor_L,lags_L] = xcorr(data1(:,1),data4(:,1));
[M_L,I_L] = max(abs(acor_L));
lagdiff14 = lags_L(I_L);
plot(lags_L, acor_L); hold on;

[acor_L,lags_L] = xcorr(data1(:,1),data5(:,1));
[M_L,I_L] = max(abs(acor_L));
lagdiff15 = lags_L(I_L);
plot(lags_L, acor_L); hold on;

data2 = circshift(data2, lagdiff12);
data3 = circshift(data3, lagdiff13);
data4 = circshift(data4, lagdiff14);
data5 = circshift(data5, lagdiff15);

ts = 1;
tf = 8;

% data1 = data1(ts*fs:tf*fs,:);
% data2 = data2(ts*fs:tf*fs,:);
% data3 = data3(ts*fs:tf*fs,:);
% data4 = data4(ts*fs:tf*fs,:);
% data5 = data5(ts*fs:tf*fs,:);

time1 = (1:length(data1))/fs;
time2 = (1:length(data2))/fs;
time3 = (1:length(data3))/fs;
time4 = (1:length(data4))/fs;
time5 = (1:length(data5))/fs;

figure(1)
plot(time1, data1(:,1))
hold on;
plot(time2, data2(:,1))
plot(time3, data3(:,1))
plot(time4, data4(:,1))
plot(time5, data5(:,1))

size(data1)
size(data2)
size(data3)
size(data4)
size(data5)

[coh_12, freq_12] = audio_mscohere(data1(:,1),data2(:,1),fs);
[coh_13, freq_13] = audio_mscohere(data1(:,1),data3(:,1),fs);
[coh_14, freq_14] = audio_mscohere(data1(:,1),data4(:,1),fs);
[coh_15, freq_15] = audio_mscohere(data1(:,1),data5(:,1),fs);
    
figure(2); hold on;
plot(freq_12, coh_12)
plot(freq_13, coh_13)
plot(freq_14, coh_14)
plot(freq_15, coh_15)
set(gca, 'XScale', 'log');
set(gca,'YGrid','on')
set(gca,'XGrid','on')
xlabel('Frequency (Hz)')
title(strcat('Coherences next groove Left Channel'))

