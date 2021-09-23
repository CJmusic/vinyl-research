function [output, info_array] = SeperateLacquer(file, offset)
                % if ismac() == true
                %     [ref, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav');
                % end
                % if ispc() == true
                %     [ref, fs] = audioread(' d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav');
                % end
        
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
                                'leadout'}% 14
                %% Reference 02072019_A0000B000r27a.wav 
                % offset = 10.625; 
                % transition = 517.375; 
                % lockoutClipped = 953.746;
                % lagdiff = []
        
                % 031418_A0000B0000r27a.wav as reference timestamps
                % offset = 10.625; 
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
                                    [383, 418]] % 12. sweep vertical 
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
         
            %~~~~~~~~~~~~~~~~~~~~~ LINE UP ~~~~~~~~~~~~~~~~~~~~~~~~% 
        
                % lockout = 950; 
                % refLockout = ref(floor(lockout*96000):end,:);
                % %% lineup audio with reference 
                % dataLockout = data(floor(950*fs):end,:);
                % disp(strcat('time diff to ref... ', num2str (length(data)/fs - length(ref)/fs)))
                % disp(strcat('size dataLockout... ', num2str(size(dataLockout))))
                % disp(strcat('size refLockout...  ', num2str(size(refLockout))))
                % fs
                % figure(2);
                % plot(refLockout);
                % hold on;
                % plot(dataLockout);
                % nblahblkah
                %% lining up audio 
                % [acor_L,lags_L] = xcorr(refLockout(:,1),dataLockout(:,1));
                % [M_L,I_L] = max(abs(acor_L));
                % lagdiff_L = lags_L(I_L);
                % lagdiff = lagdiff_L;
                lagdiff = -offset*fs;
                lagdiff = 0;
        
                % disp(strcat('lagdiff ...', num2str(lagdiff)))
        
                % % timeref = (0:length(ref)-1)/fs;
                % timedata = (0:length(data)-1)/fs  + lagdiff/fs;
                timedata = (0:length(data)-1)/fs;
        
                %***   DEBUG   ***%
                % disp(strcat(num2str(size(timeref)), num2str(size(refLockout))))
                % disp(strcat(num2str(size(timedata)), num2str(size(dataLockout))))
                
        
        
                % figure(100); grid on;
                % plot(timeref,refLockout)
                % plot(timedata,dataLockout)
                % title(track_name)
                %***   DEBUG ENDS  ***%
               
        
        %~~~~~~~~~~~~~~~~~~ NORMALIZATION RMS ~~~~~~~~~~~~~~~~~~~~%


                t = 1;
                sigtime = timedata(floor(timestamps(t,1)*fs):floor(timestamps(t,2)*fs));
                sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);
                N = 3*fs;
                seg = sig(0.33*length(sig):0.33*length(sig) + N - 1,:);


                
                disp(strcat('RMS before norm... ', num2str((rms(seg)))))
                disp(strcat('dB... ', num2str(20*log10(rms(seg)))))
                
                normalization = rms_response(seg);

                normalization_L = normalization(1);
                normalization_R = normalization(2);


                data(:,1)=data(:,1)./normalization_L;
                data(:,2)=data(:,2)./normalization_R;

                disp(strcat('normalization_L...', num2str(normalization_L)))
                disp(strcat('normalization_R...', num2str(normalization_R)))


                sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);
                N = 3*fs;
                seg = sig(0.33*length(sig):0.33*length(sig) + N - 1,:);

                disp(strcat('RMS after norm... ', num2str((rms(seg)))))
                disp(strcat('dB... ', num2str(20*log10(rms(seg)))))

    %~~~~~~~~~~~~~~~~~~ NORMALIZATION RMS ENDS ~~~~~~~~~~~~~~~~~~~~%


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
        
                        % refT = ref(1 : floor(timestamps(1,1)*fs) - lagdiff,:);
                    elseif t == length(signal_names)
                        sig = data(floor(timestamps(end,2)*fs) : length(data),:);
                        sigtime = timedata(floor(timestamps(end,2)*fs): length(data));  
        
                        % refT = ref(floor(timestamps(end,2)*fs) : length(ref),:);
                    else
                        sig = data(floor(timestamps(t-1,1)*fs): floor(timestamps(t-1,2)*fs),:);
                        sigtime = timedata(floor(timestamps(t-1,1)*fs) :floor(timestamps(t-1,2)*fs));  
        
                        % floor(timestamps(t-1,1)*fs) - lagdiff
                        % floor(timestamps(t-1,2)*fs) - lagdiff
        
                        % refT = ref(floor(timestamps(t-1,1)*fs) - lagdiff : floor(timestamps(t-1,2)*fs) - lagdiff,:);
                    end
                    % tracks(signal_names(i)) = sig;
                    signals{t} = sig;
                end
    
                output = containers.Map(signal_names, signals)
                info_array = [lagdiff, normalization_L, normalization_R];

            end