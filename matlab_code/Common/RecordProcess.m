
%~~~~~~~   TESTING   ~~~~~~~%
        
% file = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0137B0137\maxoverlap3a.wav';
% recordProcessTest(file)
%~~~~~~~ TESTING ENDS ~~~~~~~%




function output = recordProcess(file)

    % % function output = recordProcessTest(file)
    %     %~~~~~~~~~~~~~~~~~ LOAD REFERENCE ~~~~~~~~~~~~~~~~~%
    %         % try 
            % addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    %         disp('SIDE');
    %         file(length(file)-4);
    %         if ismac() == true
    %             % %~~~~ MAC ~~~~%
    %             disp('MAC')
    %             if file(length(file)-4) == 'a'
    %                 [refa, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
    %                 offseta = 15; 
    %                 [refb, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028b.wav'); 
    %                 offsetb = 13.8073; 
    %             elseif file(length(file)-4) == 'b'
    %                 [refa, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
    %                 offseta = 15;                     
    %                 [refb, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028b.wav'); 
    %                 offsetb = 13.8073; 
    %             end
    %         end 
    %         if ispc() == true
    %             %~~~~ WINDOWS ~~~~%
    %             if file(length(file)-4) == 'a'
    %                 % [ref, fs] = audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\031419_A0000B0000r028a.wav');
    %                 % offset = 10.625; 
    %                 [refa, fs] = audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\031419_A0000B0000r028a.wav');
    %                 offseta = 15; 
    %                 [refb, fs] = audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\031419_A0000B0000r028b.wav');
    %                 offsetb = 13.8073; 
    %             end
    %             if file(length(file)-4) == 'b'
    %                 % [ref, fs] = audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\031419_A0000B0000r028b.wav');
    %                 % offset = 13.8073; 
    %                 [refa, fs] = audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\031419_A0000B0000r028a.wav');
    %                 offseta = 15; 
    %                 [refb, fs] = audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\031419_A0000B0000r028b.wav');
    %                 offsetb = 13.8073; 
    %             end

    %         end
    %     %~~~~~~~~~~~~~~~~~~ Reference info ~~~~~~~~~~~~~~~~%
    
    %         timestamps_ref = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
    %         % this is how many seconds each signal is according to Chris Muths track listing
    %         lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
    %         signal_names = {'leadin',    % 1
    %                         '1kHz',      % 2
    %                         '10kHz',     % 3
    %                         '100Hz',     % 4
    %                         'sweep',     % 5
    %                         'quiet',     % 6
    %                         '3150Hz',    % 7 
    %                         '1kHzL',     % 8
    %                         'sweepL',    % 9
    %                         '1kHzR',     % 10
    %                         'sweepR',    % 11
    %                         '1kHzV',     % 12
    %                         'sweepV',    % 13
    %                         'transition',% 14
    %                         '1kHz2',     % 15
    %                         '10kHz2',    % 16
    %                         '100Hz2',    % 17
    %                         'freqsweep2',% 18
    %                         'quiet2',    % 19
    %                         '3150Hz2',   % 20
    %                         '1kHzL2',    % 21
    %                         'sweepL2',   % 22
    %                         '1kHzR2',    % 23
    %                         'sweepR2',   % 24
    %                         '1kHzV2',    % 25
    %                         'sweepV2',   % 26
    %                         'leadout'    % 27
    %         };
    %         timestamps =       [[0, 61],    % 1. 1 kHz
    %                             [61,91],    % 2. 10 kHz
    %                             [91,121],   % 3. 100 Hz
    %                             [121,159],  % 4. sweep
    %                             [159,180],  % 5. quiet
    %                             [180,245],  % 6. 3150 Hz
    %                             [245,267],  % 7. 1 kHz left
    %                             [267, 302], % 8. sweep left
    %                             [302, 325], % 9. 1 kHz right
    %                             [325, 361], % 10. sweep right
    %                             [361, 383], % 11. 1 kHz vertical
    %                             [383, 418], % 12. sweep vertical
    %                             [418, 515], % 13. transition
    %                             [515, 578], % 14. 1 kHz
    %                             [578, 608], % 15. 10 kHz
    %                             [608, 639], % 16. 100 Hz
    %                             [639, 676], % 17. sweep
    %                             [676, 698], % 18. quiet
    %                             [698, 760], % 19. 3150 Hz
    %                             [760, 785], % 20. 1 kHz left
    %                             [785, 820], % 21. sweep left
    %                             [820, 842], % 22. 1 kHz right
    %                             [842, 878], % 23. sweep right
    %                             [878, 900], % 24. 1 kHz vertical
    %                             [900, 938]];% 25. sweep vertical  
    %                             % [938, 950]];               
    %                             %% dont forget lead in and leadout
    %         timestampsa = timestamps + offseta;
    %         timestampsb = timestamps + offsetb;
    %         % timestamps = timestamps + offset;
    %         if file(length(file)-4) == 'a'
    %             timestamps = timestampsa;
    %             ref = refa;
    %         end
    %         if file(length(file)-4) == 'b'
    %             timestamps = timestampsa;
    %             ref = refb;
    %         end

    
    
    
    %     %~~~~~~~~~~~~~~~~~    LOAD FILE    ~~~~~~~~~~~~~~~~%
    
    %         [data, fs] = audioread(file);
    %         output = {};
    
    %     %~~~~~~~~~~~~~~~~~~~~~ LINE UP ~~~~~~~~~~~~~~~~~~~~~~~~% 
    
    %         lockout = 950; 
    %         refLockout = ref(floor(lockout*96000):end,:);
    %         %% lineup audio with reference 
    %         dataLockout = data(floor(950*fs):end,:);
    %         disp(strcat('time diff to ref... ', num2str (length(data)/fs - length(ref)/fs)))
    %         disp(strcat('size dataLockout... ', num2str(size(dataLockout))))
    %         disp(strcat('size refLockout...  ', num2str(size(refLockout))))
    %         fs


    %         [acor_L,lags_L] = xcorr(refLockout(:,1),dataLockout(:,1));
    %         [M_L,I_L] = max(abs(acor_L));
    %         lagdiff_L = lags_L(I_L);
    %         lagdiff = lagdiff_L;
    
    %         disp(strcat('lagdiff ...', num2str(lagdiff)))
    
    %         timeref = (0:length(ref)-1)/fs;
    %         timedata = (0:length(data)-1)/fs  + lagdiff/fs;
    
    %         %***   DEBUG   ***%
    %         % disp(strcat(num2str(size(timeref)), num2str(size(refLockout))))
    %         % disp(strcat(num2str(size(timedata)), num2str(size(dataLockout))))
            
    
    
    %         % figure(100); grid on;
    %         % plot(timeref,refLockout)
    %         % plot(timedata,dataLockout)
    %         % title(track_name)
    %         %***   DEBUG ENDS  ***%
           
    
    %     %~~~~~~~~~~~~~~~~~~~~ NORMALIZATION ~~~~~~~~~~~~~~~~~~~~%
    %         t = 1;
    %         sigtime = timedata(floor(timestamps(t,1)*fs):floor(timestamps(t,2)*fs));
    %         sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);
    
    %         % figure(t)
    %         % plot(sigtime,sig)
    
    
    %         % floor(timestamps(t,1)*fs) - lagdiff 
    %         % floor(timestamps(t,2)*fs) - lagdiff
    
    %         sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
    %         sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
    %         sigRMS=rms(sig);
    %         % normalization=sqrt(2)*sigRMS*40/7; %digital value of peak level
    %         normalization=sigRMS; %digital value of peak level
    %         data(:,1)=data(:,1)/normalization(1);% now normalized to 40cm/s peak    
    %         data(:,2)=data(:,2)/normalization(2);% now normalized to 40cm/s peak 
    %         normalization_L = normalization(1);
    %         normalization_R = normalization(2);
    
    %         disp(strcat('normalization_L... ', num2str(size(normalization_L))))
    %         disp(strcat('normalization_R... ', num2str(size(normalization_R))))
            output = [];
            fs = 96000;
            [tracks, info_array] = SeperateTracks(file);
            
            info_array
            lagdiff = info_array(1);
            normalization_L = info_array(2);
            normalization_R = info_array(3);

            signal_names = tracks.keys;
            signals = tracks.values;

            display('LOADING REFERENCE')
            if ismac() == true
                if file(length(file)-4) == 'a'
                    reference = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
                    %% Reference 02072019_A0000B000r27a.wav 
                    offset = 15; 
                elseif file(length(file)-4) == 'b'
                    disp('PC')
                    reference = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028b.wav'); 
                    %% Reference 02072019_A0000B000r27b.wav 
                    offset = 13.1;
                else 
                    disp('NO SIDE FOUND, USING SIDE A REFERENCE')
                    reference = ('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
                    %% Reference 02072019_A0000B000r27a.wav 
                    offset = 15; 
                end
            end
            if ispc() == true
                disp('IS PC')
                if file(length(file)-4) == 'a'
                    reference = ('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
                    %% Reference 02072019_A0000B000r27a.wav 
                    offset = 15; 
                elseif file(length(file)-4) == 'b'
                    disp('PC')
                    reference = ('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028b.wav'); 
                    %% Reference 02072019_A0000B000r27b.wav 
                    offset = 13.1;
                else 
                    disp('NO SIDE FOUND, USING SIDE A REFERENCE')
                    reference = ('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
                    %% Reference 02072019_A0000B000r27a.wav 
                    offset = 15; 
                end
            end
            % reference = ;
            [reftracks, ~] = SeperateTracks(reference);
            ref_signal_names = reftracks.keys;
            ref_signals = reftracks.values;

            
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
                sig = signals{t};
                refT = ref_signals{t};
                disp(strcat('track  ...',track_name))
    
    
                csig = [];
                CLICKS_R = [];
                CLICKS_L = [];
                RMS_L = [];
                RMS_R = [];
                THD_L = [];
                THD_R = [];
                
                % if t == 1
                %     sig = data(1 : floor(timestamps(1,1)*fs) - lagdiff,:);
                %     sigtime = timedata(1 : floor(timestamps(1,1)*fs) - lagdiff);  
    
                %     refT = ref(1 : floor(timestamps(1,1)*fs) - lagdiff,:);
                % elseif t == length(signal_names)
                %     sig = data(floor(timestamps(end,2)*fs) - lagdiff : length(data),:);
                %     sigtime = timedata(floor(timestamps(end,2)*fs) - lagdiff : length(data));  
    
                %     refT = ref(floor(timestamps(end,2)*fs) - lagdiff : length(ref),:);
                % else
                %     sig = data(floor(timestamps(t-1,1)*fs) - lagdiff : floor(timestamps(t-1,2)*fs) - lagdiff,:);
                %     sigtime = timedata(floor(timestamps(t-1,1)*fs) - lagdiff :floor(timestamps(t-1,2)*fs) - lagdiff);  
    
                %     % floor(timestamps(t-1,1)*fs) - lagdiff
                %     % floor(timestamps(t-1,2)*fs) - lagdiff
    
                %     refT = ref(floor(timestamps(t-1,1)*fs) - lagdiff : floor(timestamps(t-1,2)*fs) - lagdiff,:);
                % end
                
                [csig(:,1), CLICKS_L] = ClickDetect(sig(:,1));
                [csig(:,2), CLICKS_R] = ClickDetect(sig(:,2));



               
                [~, REFSa_L] = ClickDetect(refT(:,1));
                [~, REFSa_R] = ClickDetect(refT(:,2));
                [~, REFSb_L] = ClickDetect(refT(:,1));
                [~, REFSb_R] = ClickDetect(refT(:,2));
    
                % need to do the reference here by track 
                commonclicksa_L = CommonClicks(CLICKS_L, REFSa_L);
                commonclicksa_R = CommonClicks(CLICKS_R, REFSa_R);
                commonclicksb_L = CommonClicks(CLICKS_L, REFSb_L);
                commonclicksb_R = CommonClicks(CLICKS_R, REFSb_R);
    
                
                clicks_L = length(CLICKS_L);
                clicks_R = length(CLICKS_R); 

                rmssig_L = rms(sig(:,1));
                rmssig_R = rms(sig(:,2));
                rmscsig_L = rms(csig(:,1));
                rmscsig_R = rms(csig(:,2));

                RMSclicks_L = 20.0*log10(sqrt(rmssig_L^2-rmscsig_L^2));
                RMSclicks_R = 20.0*log10(sqrt(rmssig_R^2-rmscsig_R^2));

    
                RMS_L = 20.0*log10(rms(csig(:,1)));
                RMS_R = 20.0*log10(rms(csig(:,2)));
    
                Aw = audio_Aweighting(csig);
                CCIRw = audio_CCIRweighting(csig);
    
                A_L = 20.0*log10(rms_response(Aw(:,1)));
                A_R = 20.0*log10(rms_response(Aw(:,2)));
    
                CCIR_L = 20.0*log10(avg_response(CCIRw(:,1)));
                CCIR_R = 20.0*log10(avg_response(CCIRw(:,2)));
    
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
                    if strcmp(signal_names(t),'3150Hz')
                        R = 13;
                    end
                    if strcmp(signal_names(t),'3150Hz2')
                        R = 7.8;
                    end
                    centreholeoffset = (wow_L/3150)*R*cosd(20);
                else 
                    wow_L = 0;
                    wow_R = 0;
                    centreholeoffset = 0;
                end
              
                track = signal_names(t);
          
    
    
                output = [output; track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, A_L, A_R, CCIR_L, CCIR_R, clicks_L, clicks_R, commonclicksa_L, commonclicksa_R ,commonclicksb_L, commonclicksb_R, RMSclicks_L, RMSclicks_R, THD_L, THD_R, wow_L, wow_R, centreholeoffset, stereo_bleed];
               
    
            end
    end 