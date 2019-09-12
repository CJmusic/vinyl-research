close all; clear all; clc;
%% load audio
% addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000') %MAC
% [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav');

addpath('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000') %WINDOWS 
[data, fs] = audioread('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\03141_A0000B0000r30a.wav');



% figure(1)
% plot(data)

%% load reference leadout 
% [ref, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav');
[ref, fs] = audioread('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\031418_A0000B0000r27a.wav');

timestamps_ref = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
% this is how many seconds each signal is according to Chris Muth's track listing
lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
signal_names = {'leadin','1kHz', '10kHz', '100Hz', 'sweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition', '1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2', '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout'};

%% Reference 02072019_A0000B000r27a.wav 
offset = 10.625; 
transition = 517.375; 
% lockoutClipped = 953.746;
lockout = 950; 
refLockout = ref(floor(lockout*fs):end,:);
%% lineup audio with reference 
dataLockout = data(floor(950*fs):end,:);
% lagdiff = []
disp('sizes')
size(data)
size(ref)

size(refLockout)
size(dataLockout)

figure(2); 
hold on; grid on;
plot(refLockout)
plot(dataLockout)
title('pre-lineup')

[acor_L,lags_L] = xcorr(refLockout(:,1),dataLockout(:,1));
[M_L,I_L] = max(abs(acor_L));
lagdiff_L = lags_L(I_L);
lagdiff = lagdiff_L;
disp('lagdiff')
lagdiff 
datacorr = circshift(data,lagdiff);
corrLockout = datacorr(floor(lockout*fs):end,:);


figure(3); 
hold on; grid on;
plot(refLockout(:,1))
plot(corrLockout(:,1))

title('post-lineup')


figure(4); 
hold on; grid on;
plot(ref(1:20*fs,1))
plot(datacorr(1:20*fs,1))
title('post-lineup')


figure(5); 
hold on; grid on;
timeref = (0:length(ref)-1)/fs;
timedata = (0:length(data)-1)/fs  + lagdiff/fs;
plot(timeref(1:20*fs),ref(1:20*fs,1))
plot(timedata(1:20*fs),data(1:20*fs,1))
title('time-lineup')


disp('SIZES')
size(ref)
size(timeref)
size(data)
size(timedata)


figure(6); 
hold on; grid on;
plot(timeref(floor(950*fs):end),ref(floor(950*fs):end,:))
plot(timedata(floor(950*fs):end),data(950*fs:end,:))
title('pre-lineup')

% figure(4); 
% plot(data(1:20*fs,:))
% hold on;
% plot(ref(1:20*fs,:))
% grid on;
% title('post-lineup')


% get 1 kHz level and normalize 


