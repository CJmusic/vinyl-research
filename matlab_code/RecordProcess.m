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
try 
    [ref, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav');
catch
    [ref, fs] = audioread('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\031418_A0000B0000r27a.wav');
end 

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

%% PLOT WHOLE LINED UP FILES %%
% figure(1); grid on; hold on;
% plot(timeref,ref(:,1))
% plot(timedata, data(:,1))
% for xi = 1:length(timestamps)
%             x1 = timestamps(xi,1);
%             figure(1); hold on; grid on;
%             line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
% end
%% PLOTTING ENDS %%

%% PLOTTING INDIVIDUAL TRACKS %%

%~~~0. leadin 
% figure(100);
% hold on; grid on;
% plot(timedata(1:timestamps(1,1)*fs),data(1:timestamps(1,1)*fs,:));
% plot(timeref(1:timestamps(1,1)*fs),ref(1:timestamps(1,1)*fs,:));

%~~~1. 1 kHz
% get 1 kHz level
sigtime = timedata(timestamps(1,1)*fs - lagdiff :timestamps(1,2)*fs - lagdiff);
sig = data(timestamps(1,1)*fs - lagdiff : timestamps(1,2)*fs - lagdiff,:);

% normalize
sigRMS=sqrt(sum(sig)/((timestamps(1,1)-timestamps(1,2))*fs+1)); % for 7 cm/s peak
normalization=sqrt(2)*sigRMS*40/7; %digital value of peak level
data=data/normalization;% now normalized to 40cm/s peak    

% measure harmonics
THD = thd(sig(:,1),fs);
% thd(sig(:,1),fs);
% set(gca, 'XScale', 'log');


%~~~2. 10kHz 
%~~~3. 100Hz 
%~~~4. sweep 
%~~~5. quiet 
%~~~6. 3150Hz 
%~~~7. 1kHzL 
%~~~8. sweepL 
%~~~9. 1kHzR 
%~~~10. sweepR 
%~~~11. 1kHzV 
%~~~12. sweepV
%~~~13. transition 
%~~~14. 1kHz2 
%~~~15. 10kHz2 
%~~~16. 100Hz2 
%~~~17. freqsweep2 
%~~~18. quiet2 
%~~~19. 3150Hz2  
%~~~20. 1kHzL2 
%~~~21. sweepL2 
%~~~22. 1kHzR2 
%~~~23. sweepR2 
%~~~24. 1kHzV2 
%~~~25. sweepV2
%~~~26. leadout
% figure(200);
% hold on; grid on;
% plot(timedata(timestamps(25,2)*fs:end),data(timestamps(25,2)*fs:end,:));
% plot(timeref(timestamps(25,2)*fs:end),ref(timestamps(25,2)*fs:end,:));



% for i=(1:length(timestamps))
%     % if ismember(record.signal_names(i),{'1kHz','1kHzL', '1kHzR', '1kHzV','1kHz2','1kHzL2', '1kHzR2', '1kHzV2'});
%     %     disp('1kHz')
%     %     %% Flattop window to measure the peaks 
%     %     data=sqrt(sum(data)/(ntf-nts+1)); % for 7 cm/s peak
%     %     %% RMS VALUE weighted bynumber of points - CZ 
%     %     %------------signal amplitude factor for 0dB=40cm/s peak---------
%     %     normalization=sqrt(2)*rmsref*40/7; %digital value of peak level
%     %     data=data/normalization;% now normalized to 40cm/s peak    
%     % end % 1 kHz tests


%     % if ismember(record.signal_names(i),{'10kHz','10kHz', '10kHzL', '10kHzR', '10kHz2','10kHz2', '10kHzL2', '10kHzR2'});
%     %         disp('10kHz')

%     % end % 10 kHz tests


%     % if ismember(record.signal_names(i),{'100Hz', '100HzL', '100HzR','100Hz2', '100HzL2', '100HzR2'});
%     %         disp('100Hz')

%     % end % 100 Hz tests


%     % if ismember(record.signal_names(i),{'sweep', 'sweepL', 'sweepR', 'sweepV', 'sweep2', 'sweepL2', 'sweepR2', 'sweepV2'});
%     %         disp('sweep')

%     % end % sweep tests


%     % if ismember(record.signal_names(i),{'3150Hz','3150Hz2'});
%     %         disp('3150Hz')
%     % end % 3150Hz tests


%     % if ismember(record.signal_names(i),{'transition'});
%     %         disp('transition')

%     % end % transition tests


%     figure(i); hold on;
%     plot(timeref(timestamps(i,1)*fs:timestamps(i,2)*fs),ref(timestamps(i,1)*fs:timestamps(i,2)*fs,1));
%     plot(timedata(timestamps(i,1)*fs - lagdiff :timestamps(i,2)*fs - lagdiff),data(timestamps(i,1)*fs - lagdiff :timestamps(i,2)*fs - lagdiff,1));
% end

%~~~ LEADOUT ~~~%
% figure(200);
% hold on; grid on;
% plot(timedata(timestamps(25,2)*fs:end),data(timestamps(25,2)*fs:end,:));
% plot(timeref(timestamps(25,2)*fs:end),ref(timestamps(25,2)*fs:end,:));



