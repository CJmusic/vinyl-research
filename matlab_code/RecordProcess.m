close all; clear all; clc;
%% load audio
% addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000') %MAC
% [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r30a.wav');

addpath('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000') %WINDOWS 
% [data, fs] = audioread('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\03141_A0000B0000r30a.wav');

[data, fs] = audioread('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\03141_A0000B0000r28a.wav');


% figure(1)
% plot(data)

%% load reference leadout 
% [ref, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav');
[ref, fs] = audioread('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\031418_A0000B0000r27a.wav');

timestamps_ref = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
% this is how many seconds each signal is according to Chris Muth's track listing
lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
signal_names = {'leadin','1kHz', '10kHz', '100Hz', 'sweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition', '1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2',  '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout'};
%% Reference 02072019_A0000B000r27a.wav 
offset = 10.625; 
transition = 517.375; 
% lockoutClipped = 953.746;
lockout = 950; 
refLockout = ref(floor(lockout*fs):end,:);
%% lineup audio with reference 
dataLockout = data(floor(950*fs):end,:);
% lagdiff = []

% 031418_A0000B0000r27a.wav as reference timestamps
offset = 10.625; 
timestamps =       [[0, 61],   % 1 kHz
                    [61,91],    % 10 kHz
                    [91,121],   % 100 Hz
                    [121,159],  % sweep
                    [159,180],  % quiet
                    [180,245],  % 3150 Hz
                    [245,267],  % 1 kHz left
                    [267, 302], % sweep left
                    [302, 325], % 1 kHz right
                    [325, 361], % sweep right
                    [361, 383], % 1 kHz vertical
                    [383, 418], % sweep vertical
                    [418, 515], % transition
                    [515, 578], % 1 kHz
                    [578, 608], % 10 kHz
                    [608, 639], % 100 Hz
                    [639, 676], % sweep
                    [676, 698], % quiet
                    [698, 760], % 3150 Hz
                    [760, 785], % 1 kHz left
                    [785, 820], % sweep left
                    [820, 842], % 1 kHz right
                    [842, 878], % sweep right
                    [878, 900], % 1 kHz vertical
                    [900, 938]]; % sweep vertical  
                    % [938, 950]];               
                    %% dont forget lead in and leadout
timestamps = timestamps + offset;
%% lining up audio 
[acor_L,lags_L] = xcorr(refLockout(:,1),dataLockout(:,1));
[M_L,I_L] = max(abs(acor_L));
lagdiff_L = lags_L(I_L);
lagdiff = lagdiff_L;

timeref = (0:length(ref)-1)/fs;
timedata = (0:length(data)-1)/fs  + lagdiff/fs;
timestampsdata = timestamps %+ 6;

figure(1); grid on; hold on;
plot(timeref,ref(:,1))
plot(timedata, data(:,1))
for xi = 1:length(timestamps)
            x1 = timestamps(xi,1);
            figure(1); hold on; grid on;
            line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
end

%% PLOTTING INDIVIDUAL TRACKS %%
% figure(100);
% hold on; grid on;
% plot(timedata(1:timestamps(1,1)*fs),data(1:timestamps(1,1)*fs,:));
% plot(timeref(1:timestamps(1,1)*fs),ref(1:timestamps(1,1)*fs,:));

% for i=(1:length(timestamps))
%     figure(i); hold on;
%     plot(timeref(timestamps(i,1)*fs:timestamps(i,2)*fs),ref(timestamps(i,1)*fs:timestamps(i,2)*fs,1));
%     plot(timedata(timestampsdata(i,1)*fs:timestampsdata(i,2)*fs),data(timestampsdata(i,1)*fs:timestampsdata(i,2)*fs,1));
% end

% figure(200);
% hold on; grid on;
% plot(timedata(timestamps(25,2)*fs:end),data(timestamps(25,2)*fs:end,:));
% plot(timeref(timestamps(25,2)*fs:end),ref(timestamps(25,2)*fs:end,:));
%% PLOTTING INDIVIDUAL TRACKS END %%

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


