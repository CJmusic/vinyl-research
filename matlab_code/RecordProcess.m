% close all; clear all; clc;
% function [lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow, stereo_bleed] = RecordProcess(file)
function output = RecordProcess(file)
    %~~~~~~~~~~~~~~~~~ LOAD REFERENCE ~~~~~~~~~~~~~~~~~%
        try 
            [ref, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav');

        catch
            [ref, fs] = audioread('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\031418_A0000B0000r27a.wav');
        end 

    %~~~~~~~~~~~~~~~~~~ Reference info ~~~~~~~~~~~~~~~~%

        timestamps_ref = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
        % this is how many seconds each signal is according to Chris Muths track listing
        lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
        signal_names = {'leadin', '1kHz', '10kHz', '100Hz', 'sweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition', '1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2',  '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout'
        };
        %% Reference 02072019_A0000B000r27a.wav 
        offset = 10.625; 
        transition = 517.375; 
        % lockoutClipped = 953.746;
        % lagdiff = []

        % 031418_A0000B0000r27a.wav as reference timestamps
        offset = 10.625; 
        timestamps =       [[0, 61],    % 1. 1 kHz
                            [61,91],    % 2. 10 kHz
                            [91,121],   % 3. 100 Hz
                            [121,159],  % 4. sweep
                            [159,180],  % 5. quiet
                            [180,245],  % 6. 3150 Hz
                            [245,267],  % 7. 1 kHz left
                            [267, 302], % 8. sweep left
                            [302, 325], % 9. 1 kHz right
                            [325, 361], % 10. sweep right
                            [361, 383], % 11. 1 kHz vertical
                            [383, 418], % 12. sweep vertical
                            [418, 515], % 13. transition
                            [515, 578], % 14. 1 kHz
                            [578, 608], % 15. 10 kHz
                            [608, 639], % 16. 100 Hz
                            [639, 676], % 17. sweep
                            [676, 698], % 18. quiet
                            [698, 760], % 19. 3150 Hz
                            [760, 785], % 20. 1 kHz left
                            [785, 820], % 21. sweep left
                            [820, 842], % 22. 1 kHz right
                            [842, 878], % 23. sweep right
                            [878, 900], % 24. 1 kHz vertical
                            [900, 938]];% 25. sweep vertical  
                            % [938, 950]];               
                            %% dont forget lead in and leadout
        timestamps = timestamps + offset;





    %~~~~~~~~~~~~~~~~~    LOAD FILE    ~~~~~~~~~~~~~~~~%
        % try
        %     addpath(/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000) %MAC
        %     [data, fs] = audioread(/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r30a.wav);
        % catch 
        %     addpath(D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000) %WINDOWS 
        %     [data, fs] = audioread(D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\03141_A0000B0000r29a.wav);
        % end

        [data, fs] = audioread(file);
        output = {};

    %~~~~~~~~~~~~~~~~~~~~~ LINE UP ~~~~~~~~~~~~~~~~~~~~~~~~% 

        lockout = 950; 
        refLockout = ref(floor(lockout*fs):end,:);
        %% lineup audio with reference 
        dataLockout = data(floor(950*fs):end,:);

        %% lining up audio 
        [acor_L,lags_L] = xcorr(refLockout(:,1),dataLockout(:,1));
        [M_L,I_L] = max(abs(acor_L));
        lagdiff_L = lags_L(I_L);
        lagdiff = lagdiff_L;

        timeref = (0:length(ref)-1)/fs;
        timedata = (0:length(data)-1)/fs  + lagdiff/fs;

    %~~~~~~~~~~~~~~~~~~~~ NORMALIZATION ~~~~~~~~~~~~~~~~~~~~%
        t = 1;
        sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
        sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
        sigRMS=rms(sig);
        normalization=sqrt(2)*sigRMS*40/7; %digital value of peak level
        data(:,1)=data(:,1)/normalization(1);% now normalized to 40cm/s peak    
        data(:,2)=data(:,2)/normalization(2);% now normalized to 40cm/s peak 
        normalization_L = normalization(1);
        normalization_R = normalization(2);

        for t = 1:length(timestamps)
            % for each track we need: 
            %  - RMS level
            %  - Clicks
            %  - THD
            
            % for unique tracks we need 
            %  - normalization
            %  - stereo bleed (still really unsure about this test)
            %  - wow and flutter 
            track_name = signal_names{t};

            csig = [];
            CLICKS_R = [];
            CLICKS_L = [];
            RMS_L = [];
            RMS_R = [];
            THD_L = [];
            THD_R = [];

            sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
        
            sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
            
            [csig(:,1), CLICKS_L] = ClickDetect(sig(:,1));
            [csig(:,2), CLICKS_R] = ClickDetect(sig(:,2));
            
            clicks_L = length(CLICKS_L);
            clicks_R = length(CLICKS_R);

            RMS_L = 20.0*log10(rms(csig(:,1)));
            RMS_R = 20.0*log10(rms(csig(:,2)));
            THD_L = thd(csig(:,1),fs);
            THD_R = thd(csig(:,2),fs);


            if ismember(signal_names(t), {'1kHzL', 'sweepL', '1kHzL2', 'sweepL2'}) 
                ratio1 = RMS_L/RMS_R;
                %%% fft based 
                % L = 2^16;
                % seg = sig(floor(length(sig)/2) - L/2:floor(length(sig)/2) + L/2 - 1,:);
                % thdL = thd(seg(:,1))
                % thdR = thd(seg(:,2))
                % win = flattopwin(L);
                % seg = seg.*win;
        
                % fftsigL = fft(seg(:,1))/L;
                % fftsigL = fftsigL(1:L/2+1);
        
                % fftsigR = fft(seg(:,2))/L;
                % fftsigR = fftsigR(1:L/2+1);
        
                % fftfreq = fs*(0:(L/2))/L;
        
                % peakL = max(real(fftsigL));
                % peakR = max(real(fftsigR));
                % ratio2 = peakL/peakR;
                stereo_bleed = ratio1;
            elseif ismember(signal_names(t), {'1kHzR', 'sweepR', '1kHzR2', 'sweepR2'}) 
                ratio1 = RMS_R/RMS_L;
                stereo_bleed = ratio1; 
            else
                stereo_bleed = 'n/a';
            end

            if ismember(signal_names(t), {'3150Hz', '3150Hz2'})
                wow_L = WowFlutter(csig(:,1));
                wow_R = 0;%WowFlutter(csig(:,2));
            else 
                wow_L = 'n/a';
                wow_R = 'n/a';
            end
            track = signal_names(t);
            output(end+1,:) = {track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow_L, wow_R, stereo_bleed};

        end
end 