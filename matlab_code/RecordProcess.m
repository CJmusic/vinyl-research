% close all; clear all; clc;
function [lagdiff, normalization, RMS_L, RMS_R, CLICKS_L, CLICKS_R, THD_L, THD_R, wow, stereo_bleed] = RecordProcess(file)
    %%%~~~~~~~~~~~~~~~~~ LOAD REFERENCE~~~~~~~~~~~~~~~~~%%%
        try 
            [ref, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav');

        catch
            [ref, fs] = audioread('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\031418_A0000B0000r27a.wav');
        end 

    %Info about reference file

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





    %~~~~~~~~~~~~~~~~~    LOAD FILE    ~~~~~~~~~~~~~~~~%
        % try
        %     addpath(/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000) %MAC
        %     [data, fs] = audioread(/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r30a.wav);
        % catch 
        %     addpath(D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000) %WINDOWS 
        %     [data, fs] = audioread(D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\03141_A0000B0000r29a.wav);
        % end

        [data, fs] = audioread(file);


    %~~~~~~~~~~~~~~~~~    LINE UP      ~~~~~~~~~~~~~~~~% 

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

        %% PLOT WHOLE LINED UP FILES %%
        % figure(1); grid on; hold on;
        % plot(timeref,ref(:,1))
        % plot(timedata, data(:,1))
        % for xi = 1:length(timestamps)
        %             x1 = timestamps(xi,1);
        %             figure(1); hold on; grid on;
        %             line([x1 x1], get(gca, ylim),Color, black,LineStyle, --);
        % end
        %% PLOTTING ENDS %%

        %~~~~~~~~~~~~ NORMALIZATION ~~~~~~~~~~~~~~%
        t = 1;
        sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
        sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
        sigRMS=rms(sig);
        normalization=sqrt(2)*sigRMS*40/7; %digital value of peak level
        data(:,1)=data(:,1)/normalization(1);% now normalized to 40cm/s peak    
        data(:,2)=data(:,2)/normalization(2);% now normalized to 40cm/s peak 


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
                wow = 0;
            else 
                wow = 'n/a';
            end
        end
        lagdiff
        normalization
        RMS_L
        RMS_R
        clicks_L
        clicks_R
        THD_L
        THD_R
        wow
        stereo_bleed

        disp('output')
        output = {lagdiff, normalization, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow, stereo_bleed}
end %function end
%     %~~~~~~~~~~~~~~~~ 0.  leadin       ~~~~~~~~~~~~~~~~% 
%         % figure(100);
%         % hold on; grid on;
%         % plot(timedata(1:timestamps(1,1)*fs),data(1:timestamps(1,1)*fs,:));
%         % plot(timeref(1:timestamps(1,1)*fs),ref(1:timestamps(1,1)*fs,:));

%     %~~~~~~~~~~~~~~~~ 1.  1 kHz        ~~~~~~~~~~~~~~~~%

%         % get 1 kHz level
%         t = 1;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%         sigRMS=rms(sig);
%         normalization=sqrt(2)*sigRMS*40/7; %digital value of peak level
%         data(:,1)=data(:,1)/normalization(1);% now normalized to 40cm/s peak    
%         data(:,2)=data(:,2)/normalization(2);% now normalized to 40cm/s peak 

%         sig(:,1)=sig(:,1)/normalization(1);% now normalized to 40cm/s peak  
%         sig(:,2)=sig(:,2)/normalization(2);% now normalized to 40cm/s peak 

%         % measure harmonics
%         THD_L = thd(sig(:,1),fs);
%         THD_R = thd(sig(:,2),fs);

%         % click detect

%     %~~~~~~~~~~~~~~~~ 2.  10kHz        ~~~~~~~~~~~~~~~~%
%         t = 2;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%         THD_L = thd(sig(:,1),fs);
%         THD_R = thd(sig(:,2),fs);

%     %~~~~~~~~~~~~~~~~ 3.  100Hz        ~~~~~~~~~~~~~~~~%
%         t = 3;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);

%     %~~~~~~~~~~~~~~~~ 4.  sweep        ~~~~~~~~~~~~~~~~%
%         disp(sweep)
%         t = 4;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);

%         % figure(1); hold on; grid on;
%         % plotspectrum(sig, fs)
%         % set(gca, XScale, log)

%     %~~~~~~~~~~~~~~~~ 5.  quiet        ~~~~~~~~~~~~~~~~%
%         t = 5;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);

%         rmsL = 20*log10(rms(sig(:,1)));
%         rmsR = 20*log10(rms(sig(:,2)));

%     %~~~~~~~~~~~~~~~~ 6.  3150Hz       ~~~~~~~~~~~~~~~~%
%         t = 6;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 7.  1kHzL        ~~~~~~~~~~~~~~~~% 
%         disp(1kHzL)
%         t = 7;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);

%         %stereo bleed 
%         %%% rms based 
%         rmsL = rms(sig(:,1));%20*log10(rms(sig(:,1)))
%         rmsR = rms(sig(:,2));%20*log10(rms(sig(:,2)))
%         ratio1 = rmsL/rmsR

%         %%% fft based 
%         L = 2^16;
%         seg = sig(floor(length(sig)/2) - L/2:floor(length(sig)/2) + L/2 - 1,:);
%         % thdL = thd(seg(:,1))
%         % thdR = thd(seg(:,2))
%         win = flattopwin(L);
%         seg = seg.*win;

%         fftsigL = fft(seg(:,1))/L;
%         fftsigL = fftsigL(1:L/2+1);

%         fftsigR = fft(seg(:,2))/L;
%         fftsigR = fftsigR(1:L/2+1);

%         fftfreq = fs*(0:(L/2))/L;

%         peakL = max(real(fftsigL));
%         peakR = max(real(fftsigR));
%         ratio2 = peakL/peakR

%     %~~~~~~~~~~~~~~~~ 8.  sweepL       ~~~~~~~~~~~~~~~~% 
%         t = 8;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 9.  1kHzR        ~~~~~~~~~~~~~~~~%  
%         t = 9;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 10. sweepR       ~~~~~~~~~~~~~~~~%  
%         t = 10;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 11. 1kHzV        ~~~~~~~~~~~~~~~~%  
%         t = 11;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 12. sweepV       ~~~~~~~~~~~~~~~~% 
%         t = 12;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 13. transition   ~~~~~~~~~~~~~~~~%
%         t = 13;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);  

%         csig(:,1) = ClickDetect(sig(:,1));
%         csig(:,2) = ClickDetect(sig(:,2));
%         figure(1); grid on;
%         plot(sigtime,sig)
%         figure(2); grid on;
%         plot(sigtime,csig)


%     %~~~~~~~~~~~~~~~~ 14. 1kHz2        ~~~~~~~~~~~~~~~~%  
%         t = 14;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 15. 10kHz2       ~~~~~~~~~~~~~~~~%  
%         t = 15;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 16. 100Hz2       ~~~~~~~~~~~~~~~~%  
%         t = 16;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 17. freqsweep2   ~~~~~~~~~~~~~~~~%
%         t = 17;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);  
%     %~~~~~~~~~~~~~~~~ 18. quiet2       ~~~~~~~~~~~~~~~~%  
%         t = 18;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 19. 3150Hz2      ~~~~~~~~~~~~~~~~%
%         t = 19;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);   
%     %~~~~~~~~~~~~~~~~ 20. 1kHzL2       ~~~~~~~~~~~~~~~~%  
%         t = 20;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 21. sweepL2      ~~~~~~~~~~~~~~~~%  
%         t = 21;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 22. 1kHzR2       ~~~~~~~~~~~~~~~~%
%         t = 22;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);  
%     %~~~~~~~~~~~~~~~~ 23. sweepR2      ~~~~~~~~~~~~~~~~%
%         t = 23;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);  
%     %~~~~~~~~~~~~~~~~ 24. 1kHzV2       ~~~~~~~~~~~~~~~~%  
%         t = 24;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
%     %~~~~~~~~~~~~~~~~ 25. sweepV2      ~~~~~~~~~~~~~~~~%
%         t = 25;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:); 
%     %~~~~~~~~~~~~~~~~ 26. leadout      ~~~~~~~~~~~~~~~~% 
%         t = 25;
%         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
%         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);




%~~~~~~~~~~~~~~~~ AUDIO FUNCTIONS ~~~~~~~~~~~~~~~~~%
% function spec = plotspectrum(sig, fs)

%     L = length(sig);
%     % win = flattopwin(L);
%     % sig = sig.*win;

%     fftsigL = fft(sig(:,1))/L;
%     fftsigL = fftsigL(1:L/2+1);

%     fftsigR = fft(sig(:,2))/L;
%     fftsigR = fftsigR(1:L/2+1);

%     fftfreq = fs*(0:(L/2))/L;
   
%     % figure(1); grid on; hold on;
%     set(gca, XScale, log)
%     plot(fftfreq,20*log10(fftsigL));
%     plot(fftfreq,20*log10(fftsigR));
%     % figure(2); grid on; 
%     % plot(sigtime,sig)

% end







