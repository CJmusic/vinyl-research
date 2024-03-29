


% file = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r029a.wav')
% tracks = SeperateTracksTest(file)
% reference = SeperateTracksTest('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav')

% figure(1)
% plot(tracks('1kHz'))
% hold on;
% plot(reference('1kHz'))

function output = SeperateTracks(file)
% function output = SeperateTracksTest(file)
    % function output = recordProcessTest(file)
        %~~~~~~~~~~~~~~~~~ LOAD REFERENCE ~~~~~~~~~~~~~~~~~%
            % try 
            % addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    
            % %~~~~ MAC ~~~~%
            if ismac() == true
                if file(length(file)-4) == 'a'
                    [ref, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
                    %% Reference 02072019_A0000B000r27a.wav 
                    offset = 15; 
                elseif file(length(file)-4) == 'b'
                    disp('PC')
                    [ref, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028b.wav'); 
                    %% Reference 02072019_A0000B000r27b.wav 
                    offset = 13.1;
                end
            if ispc() == true
                [ref, fs] = audioread('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav');
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
                            'sweep2',% 18
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
            % offset = 10.625; 
            % transition = 517.375; 
    
            % offset = 15; %28a
            % % offset = 13.5; %28b
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
     
        %~~~~~~~~~~~~~~~~~~~~~ LINE UP ~~~~~~~~~~~~~~~~~~~~~~~~% 
    
            lockout = 950; 
            refLockout = ref(floor(lockout*96000):end,:);
            %% lineup audio with reference 
            dataLockout = data(floor(950*fs):end,:);
            disp(strcat('time diff to ref... ', num2str (length(data)/fs - length(ref)/fs)))
            disp(strcat('size dataLockout... ', num2str(size(dataLockout))))
            disp(strcat('size refLockout...  ', num2str(size(refLockout))))
            fs
            % figure(2);
            % plot(refLockout);
            % hold on;
            % plot(dataLockout);
            % nblahblkah
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
            sigtime = timedata(floor(timestamps(t,1)*fs):floor(timestamps(t,2)*fs));
            sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);
    
            % figure(t)
            % plot(sigtime,sig)
    
    
            % floor(timestamps(t,1)*fs) - lagdiff 
            % floor(timestamps(t,2)*fs) - lagdiff
    
            sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
            sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
            % sigRMS=max(sig);
            % normalization=sqrt(2)*sigRMS*40/7; %digital value of peak level
            % data(:,1)=data(:,1)/normalization(1);% now normalized to 40cm/s peak    
            % data(:,2)=data(:,2)/normalization(2);% now normalized to 40cm/s peak 
            % normalization_L = normalization(1);
            % normalization_R = normalization(2);

            L = 2^16;
        
            seg = sig(floor(length(sig)/2) - L/2:floor(length(sig)/2) + L/2 - 1,:);
        
            [b,a]=butter(2,2*100/fs,'high');% not really necessary with fft filter
            seg(:,1) = filter(b,a,seg(:,1));
            seg(:,2) = filter(b,a,seg(:,2));
        
            win = flattopwin(L);
            seg = seg.*win;
            fftsigL = fft(seg(:,1))/L;
            fftsigL = fftsigL(1:L/2+1);
            fftsigR = fft(seg(:,2))/L;
            fftsigR = fftsigR(1:L/2+1);
            fftfreq = fs*(0:(L/2))/L;
            
            peak_L = max(real(fftsigL));
            peak_R = max(real(fftsigR))

            sigRMS= [peak_L, peak_R]
            % normalization=sqrt(2)*sigRMS*40/7; %digital value of peak level
            normalization= sigRMS; %digital value of peak level
            data(:,1)=data(:,1)/normalization(1);% now normalized to 40cm/s peak    
            data(:,2)=data(:,2)/normalization(2);% now normalized to 40cm/s peak 
            normalization_L = normalization(1);
            normalization_R = normalization(2);

            signals = cell(length(signal_names),1);
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
                % tracks(signal_names(i)) = sig;
                signals{t} = sig;
            end

            output = containers.Map(signal_names, signals)

        end