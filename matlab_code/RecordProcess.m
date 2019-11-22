function output = recordProcess(file)
    %~~~~~~~~~~~~~~~~~ LOAD REFERENCE ~~~~~~~~~~~~~~~~~%
        try 
            [ref, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r28a.wav');
        catch
            [ref, fs] = audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\03141_A0000B0000r28a.wav');
        end 
    %~~~~~~~~~~~~~~~~~~ Reference info ~~~~~~~~~~~~~~~~%

        timestamps_ref = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
        % this is how many seconds each signal is according to Chris Muths track listing
        lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
        signal_names = {'leadin',    % 1
                        '1kHz',      % 2
                        '10kHz',     % 3
                        '100Hz',     % 4
                        'sweep',     % 5
                        'quiet',     % 6
                        '3150Hz',    % 7 
                        '1kHzL',     % 8
                        'sweepL',    % 9
                        '1kHzR',     % 10
                        'sweepR',    % 11
                        '1kHzV',     % 12
                        'sweepV',    % 13
                        'transition',% 14
                        '1kHz2',     % 15
                        '10kHz2',    % 16
                        '100Hz2',    % 17
                        'freqsweep2',% 18
                        'quiet2',    % 19
                        '3150Hz2',   % 20
                        '1kHzL2',    % 21
                        'sweepL2',   % 22
                        '1kHzR2',    % 23
                        'sweepR2',   % 24
                        '1kHzV2',    % 25
                        'sweepV2',   % 26
                        'leadout'    % 27
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
        disp(strcat('time diff to ref... ', num2str (length(data)/fs - length(ref)/fs)))
        disp(strcat('size dataLockout... ', num2str(size(dataLockout))))
        disp(strcat('size refLockout...  ', num2str(size(refLockout))))

        %% lining up audio 
        [acor_L,lags_L] = xcorr(refLockout(:,1),dataLockout(:,1));
        [M_L,I_L] = max(abs(acor_L));
        lagdiff_L = lags_L(I_L);
        lagdiff = lagdiff_L;

        disp(strcat('lagdiff ...', num2str(lagdiff)))

        timeref = (0:length(ref)-1)/fs;
        timedata = (0:length(data)-1)/fs  + lagdiff/fs;

        %***   DEBUG   ***%
        % disp(strcat(num2str(size(timeref)), num2str(size(refLockout))))
        % disp(strcat(num2str(size(timedata)), num2str(size(dataLockout))))
        


        % figure(100); grid on;
        % plot(timeref,refLockout)
        % plot(timedata,dataLockout)
        % title(track_name)
        %***   DEBUG ENDS  ***%


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

        for t = 1:length(signal_names)
            % for each track we need: 
            %  - RMS level
            %  - Clicks
            %  - THD
            
            % for unique tracks we need 
            %  - normalization
            %  - stereo bleed (still really unsure about this test)
            %  - wow and flutter 
            track_name = signal_names{t};
            disp(strcat('track  ...',track_name))


            csig = [];
            CLICKS_R = [];
            CLICKS_L = [];
            RMS_L = [];
            RMS_R = [];
            THD_L = [];
            THD_R = [];
            
            if t == 1
                sig = data(1 : floor(timestamps(1,1)*fs) - lagdiff,:);
                sigtime = timedata(1 : floor(timestamps(1,1)*fs) - lagdiff);  

                refT = ref(1 : floor(timestamps(1,1)*fs) - lagdiff,:);
            elseif t == length(signal_names)
                sig = data(floor(timestamps(end,2)*fs) - lagdiff : length(data),:);
                sigtime = timedata(floor(timestamps(end,2)*fs) - lagdiff : length(data));  

                refT = ref(floor(timestamps(end,2)*fs) - lagdiff : length(ref),:);
            else
                sig = data(floor(timestamps(t-1,1)*fs) - lagdiff : floor(timestamps(t-1,2)*fs) - lagdiff,:);
                sigtime = timedata(floor(timestamps(t-1,1)*fs) - lagdiff :floor(timestamps(t-1,2)*fs) - lagdiff);  

                % floor(timestamps(t-1,1)*fs) - lagdiff
                % floor(timestamps(t-1,2)*fs) - lagdiff

                refT = ref(floor(timestamps(t-1,1)*fs) - lagdiff : floor(timestamps(t-1,2)*fs) - lagdiff,:);
            end
            
            [csig(:,1), CLICKS_L] = ClickDetect(sig(:,1),200,20);
            [csig(:,2), CLICKS_R] = ClickDetect(sig(:,2),200,20);
           
            [~, REFS_L] = ClickDetect(ref(:,1),200,20);
            [~, REFS_R] = ClickDetect(ref(:,2),200,20);

            % need to do the reference here by track 
            [diff_arrayL, ~] = audio_clickmatrix(CLICKS_L, REFS_L);
            [diff_arrayR, ~] = audio_clickmatrix(CLICKS_R, REFS_R);

            commonclicks_L = length(diff_arrayL);
            commonclicks_R = length(diff_arrayL);
            
            clicks_L = length(CLICKS_L);
            clicks_R = length(CLICKS_R); 

            RMS_L = 20.0*log10(rms(csig(:,1)));
            RMS_R = 20.0*log10(rms(csig(:,2)));

            %***    DEBUG    ***%
            % figure(t); grid on;
            % plot(csig)
            % title(track_name)
            %*** DEBUG ENDS ***%


            THD_L = thd(csig(:,1),fs);
            THD_R = thd(csig(:,2),fs);


            if ismember(signal_names(t), {'1kHzL', 'sweepL', '1kHzL2', 'sweepL2'})             
                stereo_bleed = StereoBleed(csig,1);
            elseif ismember(signal_names(t), {'1kHzR', 'sweepR', '1kHzR2', 'sweepR2'}) 
                stereo_bleed = StereoBleed(csig,2);
            else
                stereo_bleed = 0;
            end

            if ismember(signal_names(t), {'3150Hz', '3150Hz2'})
                wow_L = WowFlutter(csig(:,1));
                wow_R = WowFlutter(csig(:,2));
            else 
                wow_L = 0;
                wow_R = 0;
            end
            % figure(1000+t); hold on; grid on;
            % title(signal_names(t))
            % plot(csig(:,1))
            % plot(csig(:,2))
            track = signal_names(t);
            lagdiff =num2str(lagdiff);
            normalization_L =num2str(normalization_L);
            normalization_R =num2str(normalization_R);
            RMS_L =num2str(RMS_L);
            RMS_R =num2str(RMS_R);
            clicks_L =num2str(clicks_L);
            clicks_R =num2str(clicks_R);
            THD_L =num2str(THD_L);
            THD_R =num2str(THD_R);
            wow_L =num2str(wow_L);
            wow_R =num2str(wow_R);
            stereo_bleed = num2str(stereo_bleed);

            output = [track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, commonclicks_L, commonclicks_R  THD_L, THD_R, wow_L, wow_R, stereo_bleed];
            % class(output)
            % disp('END RECORD PROCESS')

        end
end 