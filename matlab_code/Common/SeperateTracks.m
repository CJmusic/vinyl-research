% clear all;close all;clc
% set(0,'DefaultLineLineWidth',1.0);
% set(0,'DefaultAxesFontSize',12);
% set(0,'DefaultAxesFontWeight','bold')
% set(0,'DefaultAxesLineWidth',1.5)

% addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')

% tracks = SeperateTracksTest('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r029a.wav')
% reference = SeperateTracksTest('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav')

% figure(1)
% plot(tracks('1kHz'))
% hold on;
% title('normalized signal')


% fs = 96000;
% seg = tracks('1kHz');
% % seg = seg(1:2^16,:);
% L = 2^16;
% seg = seg(floor(length(seg)/2) - L/2:floor(length(seg)/2) + L/2 - 1,:);
% size(seg)

% refseg = reference('1kHz');
% refseg = refseg(1:2^16,:);


% [spec, fftfreq] = audio_spectrum(seg, fs, 1, 2^16);
% [refspec, fftfreq] = audio_spectrum(refseg, fs, 1, 2^16);

% figure(2)
% audio_plotspectrum(fftfreq, spec, 'after seperate tracks')
% % hold on;
% figure(3)
% audio_plotspectrum(fftfreq, refspec, 'reference')




function [output, info_array] = SeperateTracks(file)
% function [output, info_array] = SeperateTracksTest(file)
        %~~~~~~~~~~~~~~~~~ LOAD REFERENCE ~~~~~~~~~~~~~~~~~%
            % try 
            addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
            disp('SEPERATE TRACKS CALLED')
           
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
                            'sweep2',    % 18
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
            % timestamps =       [[0, 61],    % 1. 1 kHz
            %                     [61,91],    % 2. 10 kHz
            %                     [91,121],   % 3. 100 Hz
            %                     [121,159],  % 4. sweep
            %                     [159,180],  % 5. quiet
            %                     [180,245],  % 6. 3150 Hz
            %                     [245,267],  % 7. 1 kHz left
            %                     [267, 302], % 8. sweep left
            %                     [302, 325], % 9. 1 kHz right
            %                     [325, 361], % 10. sweep right
            %                     [361, 383], % 11. 1 kHz vertical
            %                     [383, 418], % 12. sweep vertical
            %                     [418, 515], % 13. transition
            %                     [515, 578], % 14. 1 kHz
            %                     [578, 608], % 15. 10 kHz
            %                     [608, 639], % 16. 100 Hz
            %                     [639, 676], % 17. sweep
            %                     [676, 698], % 18. quiet
            %                     [698, 760], % 19. 3150 Hz
            %                     [760, 785], % 20. 1 kHz left
            %                     [785, 820], % 21. sweep left
            %                     [820, 842], % 22. 1 kHz right
            %                     [842, 878], % 23. sweep right
            %                     [878, 900], % 24. 1 kHz vertical
            %                     [900, 938]];% 25. sweep vertical  
            %                     % [938, 950]];               
            %                     %% dont forget lead in and leadout

        timestamps =            [[14.624, 75.229],    % 1. 1 kHz
                                [75.229, 104.512],    % 2. 10 kHz
                                [104.512, 135.275],   % 3. 100 Hz
                                [135.275, 173.224],  % 4. sweep
                                [173.224, 195.111],  % 5. quiet
                                [195.111, 258.293],  % 6. 3150 Hz
                                [258.293, 281.625],  % 7. 1 kHz left
                                [281.625, 318.352], % 8. sweep left
                                [318.352, 339.302], % 9. 1 kHz right
                                [339.302, 376.067], % 10. sweep right
                                [376.067, 396.919], % 11. 1 kHz vertical
                                [396.919, 432.435], % 12. sweep vertical
                                [432.435, 532.324], % 13. transition
                                [532.324, 592.655], % 14. 1 kHz
                                [592.655, 623.211], % 15. 10 kHz
                                [623.211, 653.123], % 16. 100 Hz
                                [653.123, 690.572], % 17. sweep
                                [690.572, 712.537], % 18. quiet
                                [712.537, 775.368], % 19. 3150 Hz
                                [775.368, 799.286], % 20. 1 kHz left
                                [799.286, 835.582], % 21. sweep left
                                [835.582, 857.235], % 22. 1 kHz right
                                [857.235, 893.610], % 23. sweep right
                                [893.610, 914.697], % 24. 1 kHz vertical
                                [914.697, 950.037]];% 25. sweep vertical  
                                % [938, 950]];               
                                %% dont forget lead in and leadout


    
            [data, fs] = audioread(file);
     
        %~~~~~~~~~~~~~~~~~~~~~ LINE UP ~~~~~~~~~~~~~~~~~~~~~~~~% 
            %~~~~~~~~~~ CORRELATION ~~~~~~~~~%
            % lockout = 950; 
            % refLockout = ref(floor(lockout*96000):end,:);
            % %% lineup audio with reference 
            % dataLockout = data(floor(950*fs):end,:);
            % disp(strcat('time diff to ref... ', num2str(length(data)  - length(ref))))
            % disp(strcat('size dataLockout... ', num2str(size(dataLockout))))
            % disp(strcat('size refLockout...  ', num2str(size(refLockout))))

            % [acor_L,lags_L] = xcorr(refLockout(:,1),dataLockout(:,1));
            % [M_L,I_L] = max(abs(acor_L));
            % lagdiff_L = lags_L(I_L);
            % lagdiff = lagdiff_L;
            % disp(strcat('lagdiff...', num2str(lagdiff)))

            % timeref = (0:length(ref)-1) ;
            % timedata = (0:length(data)-1)   + lagdiff ;
            %~~~~~~~~~~ CORRELATION  ENDS ~~~~~~~~~%
    

            %~~~~~~~~~~ MANUAL LINEUP ~~~~~~~~~~~%
            timestring = file(end-11:end-4); % get the 

            timepip = str2num(timestring(1:2))*60 + str2num(timestring(3:end));
            
            timestringref = '1558.066';
            timepipref = str2num(timestringref(1:2))*60 + str2num(timestringref(3:end));



            timediff = timepipref - timepip; 

            lagdiff = floor(timediff*96000);

            %~~~~~ Correlation correction ~~~~~%

            ref = audioread('/Volumes/AUDIOBANK/audio_files/A0000B0000/031418_A0000B0000r028a1558.066/leadout.wav');
            timepipref2 = 7.934;

            lockout = timepipref2; 
            refLockout = ref(timepipref2*fs-0.25*fs:timepipref2*fs+0.25,:);
            dataLockout = data(timepip*fs-0.25*fs:timepip*fs+0.25*fs,:);

            [acor_L,lags_L2] = xcorr(refLockout(:,1),dataLockout(:,1));
            [M_L,I_L] = max(abs(acor_L));
            lagdiff_L2 = lags_L2(I_L);
            lagdiff2 = lagdiff_L2;
            lagdiff2 = 0;
            timediff = timediff + lagdiff2/fs;
            %~~~~~ Correlation correction ends ~~~~~%

            timedata = (0:length(data)-1);% + timediff;
            timestamps = timestamps - timediff;
            if timestamps(1,1) < 0; 
                timestamps(1,1) = 1;
            end

            disp(strcat('timediff...',num2str(timediff)))
            %~~~~~~~~~~ MANUAL LINEUP  ENDS ~~~~~~~~~~~%

            %~~~~~~~~~~ DEBUG ~~~~~~~~~~%
            % figure(1)
            % plot(ref(timepipref2*fs-0.25*fs:timepipref2*fs+0.25*fs))
            % hold on; grid on;
            % plot(data(timepip*fs-0.25*fs:timepip*fs+0.25*fs))



            %~~~~~~~~~~ DEBUG ENDS ~~~~~~~~~~%
            


        %~~~~~~~~~~~~~~~~~~~~ NORMALIZATION ~~~~~~~~~~~~~~~~~~~~%

            %~~~~ separate out the 1 kHz track
            t = 1;
            sigtime = timedata(floor(timestamps(t,1)*fs):floor(timestamps(t,2)*fs));
            sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);

            N = 3*fs;
            win = flattopwin(N);
            windowfactor = 0.2155774;
            
            seg = sig(0.33*length(sig):0.33*length(sig) + N - 1,:);
            winseg = seg.*win/windowfactor;
            
            [seg_fft, freq_fft] = audio_spectrum(winseg, fs, 1, N); 
            normalization = max(abs(seg_fft))/(sqrt(2));
            segnorm = seg./normalization;
            winsegnorm = segnorm.*win/windowfactor;
            [seg_fftnorm, freq_fft] = audio_spectrum(winsegnorm, fs, 1, N);
            
            normalization_L = normalization(1);
            normalization_R = normalization(2);



            data(:,1)=data(:,1)./normalization_L;
            data(:,2)=data(:,2)./normalization_R;

            disp(strcat('normalization_L...', num2str(normalization_L)))
            disp(strcat('normalization_R...', num2str(normalization_R)))
            disp(strcat('max fft before norm... ', num2str(max(abs(seg_fft)))))
            disp(strcat('dB... ', num2str(20*log10(max(abs(seg_fft))))))
            disp(strcat('max fft after norm... ', num2str(max(abs(seg_fftnorm)))))
            disp(strcat('dB... ', num2str(20*log10(max(abs(seg_fftnorm))))))
            

            sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);
            N = 3*fs;
            win = flattopwin(N);
            windowfactor = 0.2155774;
            seg = sig(0.33*length(sig):0.33*length(sig) + N - 1,:);
            winseg = seg.*win/windowfactor;
            [seg_fft, freq_fft] = audio_spectrum(winseg, fs, 1, N); 
            disp(strcat('data amplitude after norm... ', num2str(max(abs(seg_fft)))))
            disp(strcat('dB... ', num2str(20*log10(max(abs(seg_fft))))))


            figure(1)
            plot(seg)

            %~~~~~~~~~~~~~~~~~~~~~~~~~~~ END NORMALIZATION~~~~~~~~~~~~~~~~~~~~~~~~~~%



            % normalization_L = 0;
            % normalization_R = 0;

            signals = cell(length(signal_names),1);
            signal_times = cell(length(signal_names),1);

            for t = (1:length(signal_names))
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
                    sig = data(1 : floor(timestamps(1,1)*fs),:);
                    sigtime = timedata(1 : floor(timestamps(1,1)*fs)); 
                    % refT = ref(1 : floor(timestamps(1,1)*fs),:);
                elseif t == length(signal_names)
                    sig = data(floor(timestamps(end,2)*fs) : length(data),:);
                    sigtime = timedata(floor(timestamps(end,2)*fs) : length(data));  
    
                    % refT = ref(floor(timestamps(end,2)*fs) : length(ref),:);
                else
                    sig = data(floor(timestamps(t-1,1)*fs) : floor(timestamps(t-1,2)*fs),:);
                    sigtime = timedata(floor(timestamps(t-1,1)*fs) :floor(timestamps(t-1,2)*fs));  
    
                    % floor(timestamps(t-1,1) )
                    % floor(timestamps(t-1,2) )
                    % refT = ref(floor(timestamps(t-1,1)*fs) : floor(timestamps(t-1,2)*fs),:);
                end

            % tracks(signal_names(i)) = sig;
                signals{t} = sig;
                signal_times{t} = sigtime; % not currently assigned to output

                % figure(500+t)
                % plot(signals{t})
                % for xi = 1:length(clicks)
                %     x1 = (clicks(xi));
                %     figure(1); hold on;
                %     line([x1 x1], get(gca, 'ylim'),'Color', 'red','LineStyle', '--');
                % end
            
            end

            disp('ASSIGNING OUTPUT')
            output = containers.Map(signal_names, signals);
            info_array = [lagdiff, normalization_L, normalization_R];
            disp('EXITING SEPERATE TRACKS')


            %~~~~~~~~~~ DEBUG ~~~~~~~~~~%

            % timelockout = sigtime
            % timelockoutref = 

            % sig = data(floor(timestamps(end,2)*fs)  : length(data),:);
            % sigtime = timedata(floor(timestamps(end,2)*fs)  : length(data));  

            % refT = ref(floor(timestamps(end,2)*fs) : length(ref),:);


            % figure(100)
            % sigtime = (1:length(sig(:,1)))/fs;
            % plot(sigtime, sig(:,1))
            % hold on; grid on;
            % reftime = (1:length(refT(:,1)))/fs;
            % plot(reftime,refT(:,1))

            %~~~~~~~~~~ DEBUG ENDS ~~~~~~~~~~%



            % figure(10000)
            % plot(signals{2})
            % figure(10001)
            % plot(output('1kHz'))

            % figure(1001)
            % timeref  = (0:length(ref)-1)/fs;
            % timedata = (0:length(data)-1)/fs + lagdiff/fs;
            % subplot(2,1,1)
            % plot(timeref, ref(:,1),'k')
            % grid on;
            % title('Reference track')
            % xlabel('time [s]')
            % subplot(2,1,2)
            % plot(timedata,data(:,1),'k')
            % grid on;
            % title('Recorded track')
            % xlabel('time [s]')

        end