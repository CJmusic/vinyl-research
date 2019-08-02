%% This file will measure the coherences in grooves of a record 



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


filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r27a-mac.wav'; RECORD = true; 

if RECORD == true; 
    referencePath = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r28a.wav' 
    recordPath = filename; 

    reference = audio_recordclass(referencePath);
    reference.process_tracks()
    record = audio_recordclass(recordPath);
    record.process_tracks()

    record.lagdiff = -1*audio_corrlineup(record.tracks('leadout'), reference.tracks('leadout'));
    record.lagcorrect()

    data = record.tracks('transition');
    fs = record.fs;

else 
    nstart = 0;
    nend = 0; 
    tstart = 0;
    tend = 0;
    data = 0; 
end 


grLen = 1.8*fs; 
nGrooves = floor(length(data)/grLen);
grArray = [];
grTime = (0:grLen-1)/fs;

for ng = 1:nGrooves
    grArray(:,:,ng) = data(1+(ng-1)*grLen:ng*grLen,:);
end

grCorr = [];

for ng = 1:nGrooves
    % ORIGINAL WAVEFORM 
    figure(1); hold on;
    % plot(grTime,grArray(:,1,ng), 'DisplayName',['groove',num2str(ng)]) 
    plot(grTime,grArray(:,1,ng), 'DisplayName',['groove',num2str(ng)], 'Color', [1 - 1.0*ng/100 ,0,1.0*ng/100,0.5]);
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('time (s)')
    ylabel('Amplitude')

    % ORIGINAL FFT 


    % CORRELATION TO FIRST GROOVE 
    grCorr = [grCorr, sum(grArray(:,:,1)*grArray(:,:,ng))];

    % COHERENCE TO FIRST GROOVE 


    % COHERENCE TO NEXT GROOVE 


    % XCORR for lineup


end





















% %%%~~~~~~~~~~~~~ OLD CODE ~~~~~~~~~~~~~%%%


% % [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/030419_r26fivetrials/one.wav');
% % [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/042619_lacquerApostsweepv.wav');
% % [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/lacquer+pioneer/052319_pioneer.wav');
% % [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/cut_silence/031418_A0000B0000r29a-cutsilence.wav');
% [data, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r27a-mac.wav'); tstart = 435; tend = 455;%535;



% %signal_names = {'leadin','1kHz', '10kHz', '100Hz', 'freqsweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition', '1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2', '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout'};
% %timestamps = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
% %lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz

% % 1 khz
% % tstart = 3;
% % tend = 10;
% data = data(tstart*fs:tend*fs,:);

% % quiet track on record
% %tstart = 180;
% %tend = 246;

% time = linspace(0,(length(data)-1)/fs,length(data));
% rotation_speed = 33.33333;%45;
% % T = 60/rotation_speed; %this is the length of one groove groove
% grLen = 172800;
% % num_segs = (floor(length(data)/fs/T))
% num_segs = (floor(length(data)/grLen))
% % grLen = round(T*fs)
% grTime = time(1:grLen);
% grArray = []; %need to 


% for ng = 1:num_segs
%     grArray(:,:,ng) = data(1+(ng-1)*grLen:ng*grLen,:);
% end
% for ng = 1:num_segs-1

%     %%% coherences in the left channel 
%     disp('NEW LOOP~~~~~~~~~~~~~~~~~~~~~~~~~')
%     %%********* PRINT THE LAGDIFF BETWEEN GROOVES ************
%     if ng > 1;
%         groovediff = audio_corrlineup(grArray(:,1,ng), grArray(:,1,ng-1))
%         grArray(:,2,ng-1)= circshift(grArray(:,2,ng-1), groovediff);
%     end

%     %%%%%
%     [coh_nextL, freq_coh] = audio_mscohere(grArray(:,1,ng), grArray(:,1,ng+1), fs);
%     [coh_firstL, ~] = audio_mscohere(grArray(:,1,1), grArray(:,1,ng+1), fs);

%     %%% coherences in the right channel    
%     [coh_nextR, ~] = audio_mscohere(grArray(:,2,ng), grArray(:,2,ng+1), fs);
%     [coh_firstR, ~] = audio_mscohere(grArray(:,2,1), grArray(:,2,ng+1), fs);


%     figure(1); hold on;
%     plot(grTime,grArray(:,1,ng), 'DisplayName',['groove',num2str(ng)]) 
%     % plot(grTime,grArray(:,1,ng), 'DisplayName',['groove',num2str(ng)], 'Color', [1 - 1.0*ng/100 ,0,1.0*ng/100,0.5]);
%     set(gca,'YGrid','on')
%     set(gca,'XGrid','on')
%     xlabel('time (s)')
%     ylabel('Amplitude')

%     figure(2); hold on;
%     plot(freq_coh, coh_nextL, 'DisplayName',['groove',num2str(ng)])
%     set(gca, 'XScale', 'log');
%     set(gca,'YGrid','on')
%     set(gca,'XGrid','on')
%     xlabel('Frequency (Hz)')
%     title(strcat('Coherences next groove Left Channel'))


%     figure(3); hold on;
%     plot(freq_coh, coh_nextR, 'DisplayName',['groove',num2str(ng)])
%     set(gca, 'XScale', 'log');
%     set(gca,'YGrid','on')
%     set(gca,'XGrid','on')
%     xlabel('Frequency (Hz)')
%     title(strcat('Coherences next groove Right Channel'))


%     figure(4); hold on;
%     plot(freq_coh, coh_firstL, 'DisplayName',['groove',num2str(ng)])
%     set(gca, 'XScale', 'log');
%     set(gca,'YGrid','on')
%     set(gca,'XGrid','on')
%     xlabel('Frequency (Hz)')
%     title(strcat('Coherences first groove Left Channel'))

%     figure(5); hold on;
%     plot(freq_coh, coh_firstR, 'DisplayName',['groove',num2str(ng)])
%     set(gca, 'XScale', 'log');
%     set(gca,'YGrid','on')
%     set(gca,'XGrid','on')
%     xlabel('Frequency (Hz)')
%     title(strcat('Coherences first groove Right Channel'))


%     figure(6); hold on; %% plot the spectrum of each groove
%     grLen = length(grArray(:,1,ng));
%     freq_fft = fs*(0:(grLen/2))/grLen;
%     data_fft = fft(grArray(:,1,ng))/grLen;
%     data_fft = data_fft(1:grLen/2+1);

%     plot(freq_fft, 20.0*log10(data_fft)) 
%     grid on; 
%     set(gca, 'XScale', 'log');
%     title(strcat("Groove's Spectrum"))
%     xlabel('Frequency (Hz)')
%     ylabel('Level (dB)')  

%     [cor_first, lags_first] = xcorr(grArray(:,1,1), grArray(:,1,ng));
%     [cor_next, lags_next] = xcorr(grArray(:,1,ng), grArray(:,1,ng+1));
    
%     figure(7); hold on;
%     plot(lags_first, cor_first);
%     grid on; 
%     title('xcorr with first groove')
%     xlabel('lag')
%     ylabel('correlation')

%     figure(8); hold on;
%     plot(lags_next, cor_next);
%     grid on; 
%     title('xcorr with next groove')
%     xlabel('lag')
%     ylabel('correlation')
% end

% figure(1)
% legend()

% figure(2)
% legend()

% figure(3)
% legend()

% figure(4)
% legend()

% figure(5)
% legend()

% figure(6)
% legend()

% figure(7)
% legend()

% figure(8)
% legend()