close all; clear all; clc;
%% load audio
% addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000') %MAC
% [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r30a.wav');

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

timestamps =       [[0, 62],   % 1 kHz
                    [62,92],    % 10 kHz
                    [92,124],   % 100 Hz
                    [124,160],  % sweep
                    [160,182],  % quiet
                    [182,248],  % 3150 Hz
                    [248,268],  % 1 kHz left
                    [268, 306], % sweep left
                    [306, 326], % 1 kHz right
                    [326, 364], % sweep right
                    [364, 384], % 1 kHz vertical
                    [384, 421], % sweep vertical
                    [421, 517], % transition
                    [517, 579], % 1 kHz
                    [579, 609], % 10 kHz
                    [609, 641], % 100 Hz
                    [641, 677], % sweep
                    [677, 699], % quiet
                    [699, 765], % 3150 Hz
                    [765, 785], % 1 kHz left
                    [785, 823], % sweep left
                    [823, 843], % 1 kHz right
                    [843, 881], % sweep right
                    [881, 901], % 1 kHz vertical
                    [901, 938]] % sweep vertical
                    %% dont forget lead in and leadout


%% lining up audio 
[acor_L,lags_L] = xcorr(refLockout(:,1),dataLockout(:,1));
[M_L,I_L] = max(abs(acor_L));
lagdiff_L = lags_L(I_L);
lagdiff = lagdiff_L;
% datacorr = circshift(data,lagdiff);
% corrLockout = datacorr(floor(lockout*fs):end,:);
timeref = (0:length(ref)-1)/fs;
timedata = (0:length(data)-1)/fs  + lagdiff/fs;

disp('testing timestamps')
timestamps(1,1)
figure(1);
hold on; grid on;
y = data(timestamps(1,1)+1:timestamps(1,2)*fs);
x = timedata(timestamps(1,1)+1,timestamps(1,2)*fs);
plot(x,y)
% plot(timedata(timestamps(1,1)*fs:timestamps(1,2)*fs), 
%      data(timestamps(1,1)*fs:timestamps(1,2)*fs))
% plot(timeref(timestamps(1,1):timestamps(1,2)*fs), 
%      ref(timestamps(1,1)*fs:timestamps(1,2)*fs))



%%% LINEUP TEST PLOTS %%%
% figure(2); 
% hold on; grid on;
% plot(refLockout)
% plot(dataLockout)
% title('pre-lineup')


% figure(3); 
% hold on; grid on;
% plot(refLockout(:,1))
% plot(corrLockout(:,1))

% title('post-lineup')


% figure(4); 
% hold on; grid on;
% plot(ref(1:20*fs,1))
% plot(datacorr(1:20*fs,1))
% title('post-lineup')


% figure(5); 
% hold on; grid on;
% timeref = (0:length(ref)-1)/fs;
% timedata = (0:length(data)-1)/fs  + lagdiff/fs;
% plot(timeref(1:20*fs),ref(1:20*fs,1))
% plot(timedata(1:20*fs),data(1:20*fs,1))
% title('time-lineup')

% figure(6); 
% hold on; grid on;
% plot(timeref(floor(950*fs):end),ref(floor(950*fs):end,1))
% plot(timedata(floor(950*fs):end),data(950*fs:end,1))
% title('time-lineup')

% figure(7); 
% hold on; grid on;
% plot(timeref(floor(670*fs):floor(700*fs)),ref(floor(670*fs):floor(700*fs),1))
% plot(timedata(floor(670*fs):floor(700*fs)),data(floor(670*fs):floor(700*fs),1))
% title('time-lineup')
%%% LINEUP TEST PLOTS ENDS %%%


% figure(4); 
% plot(data(1:20*fs,:))
% hold on;
% plot(ref(1:20*fs,:))
% grid on;
% title('post-lineup')


% get 1 kHz level and normalize 


