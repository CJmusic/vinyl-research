clear all;close all;clc
set(0,'DefaultLineLineWidth',0.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/Code/vinyl-research/matlab_code/from_John')

timestring = '161.695'

% tracks = SeperateTracksTest('/Volumes/AUDIOBANK/audio_files/A0137B0137/039b.wav', timestring)
tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/039b1601.695.wav')

signal_names = tracks.keys;
signals = tracks.values;
fs = 96000;
for i = 1:length(signal_names);
    trackname = signal_names{i};
    sig = signals{i}; 
    time = (1:length(sig))/fs;
    plot_track(time,sig,trackname,trackname)
    
end



function plot_track(x, y, trackname, filename)
  
    fig = figure('Visible', 'off')
    subplot(2,1,1)
    title(strcat(trackname, ' left channel'))
    plot(x,y(:,1),'k')
    grid on;
    subplot(2,1,2)
    plot(x,y(:,2),'k')
    title(strcat(trackname, ' right channel'))
    grid on;
    ylim([-1,1])
    xlabel('time [s]')
    ylabel('signal level')
    plotname = strcat('Plots/ManualLineup/', filename, '.png')
    saveas(fig, plotname);




end

function [output, info_array] = SeperateTracksTest(file, timestring)
    % function [output, info_array] = SeperateTracksTest(file)
            %~~~~~~~~~~~~~~~~~ LOAD REFERENCE ~~~~~~~~~~~~~~~~~%
                % try 
                % addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
        
                % % %~~~~ MAC ~~~~%
                % disp('SEPERATE TRACKS CALLED')
                % % USE 003a and 052b as references 
    
                % if ismac() == true
                %     if file(length(file)-4) == 'a'
                %         disp('MAC')
                %         disp('Using a side reference')
                %         % [ref, ] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); offset = 15; \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ \\\\\\\\\\\\ 
                %         [ref, ] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/040318_A0000B0000r003a.wav'); offset = 9.227; 
                %     elseif file(length(file)-4) == 'b'
                %         disp('MAC')
                %         disp('Using b side reference')
    
                %         % [ref, ] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r029b.wav'); offset = 13.1;
                %         [ref, ] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/080619-A0000B0000r052b.wav'); offset = 14.662;
                %     else 
                %         disp('NO SIDE FOUND, USING SIDE A REFERENCE')
                    [ref, ] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); offset = 14.5; refpip = 1558.006;
                %         [ref, ] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/040318_A0000B0000r003a.wav'); offset = 9.227; 
                        
                %     end
                % end
                % if ispc() == true
                %     disp('IS PC')
                %     if file(length(file)-4) == 'a'
                %         % [ref, ] = audioread('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); offset = 15; 
                %         [ref, ] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/040318_A0000B0000r003a.wav'); offset = 9.227; 
    
                %     elseif file(length(file)-4) == 'b'
                %         disp('PC')
                %         % [ref, ] = audioread('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028b.wav'); offset = 13.1;
                %         [ref, ] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/080619-A0000B0000r052b.wav'); offset = 14.662;
                %     else 
                %         disp('NO SIDE FOUND, USING SIDE A REFERENCE')
                %         % [ref, ] = audioread('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); offset = 15; 
                %         [ref, ] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/040318_A0000B0000r003a.wav'); offset = 9.227; 
                %     end
    
                % end
                
        
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
                % offset = 10.625; 
                % transition = 517.375; 
        
                % offset = 15; %28a
                % % offset = 13.5; %28b
                timestamps =       [[0, 61],    % 1. 1 kHz
                                    [61,90],    % 2. 10 kHz
                                    [90,121],   % 3. 100 Hz
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
        
        
        
        
        
                [data, fs] = audioread(file);
         
            %~~~~~~~~~~~~~~~~~~~~~ LINE UP ~~~~~~~~~~~~~~~~~~~~~~~~% 
        
                % lockout = 950; 
                % refLockout = ref(floor(lockout*96000):end,:);
                % %% lineup audio                                                     with reference 
                % dataLockout = data(floor(950*fs):end,:);
                % disp(strcat('time diff to ref... ', num2str(length(data)  - length(ref))))
                % disp(strcat('size dataLockout... ', num2str(size(dataLockout))))
                % disp(strcat('size refLockout...  ', num2str(size(refLockout))))
    
                % [acor_L,lags_L] = xcorr(refLockout(:,1),dataLockout(:,1));
                % [M_L,I_L] = max(abs(acor_L));
                % lagdiff_L = lags_L(I_L);
                % lagdiff = lagdiff_L;
                % % figure(21)
                % % plot(lags_L/fs, acor_L)
    
                % disp(strcat('lagdiff...', num2str(lagdiff)))
                timediff = str2num(timestring(1:2))*60 + str2num(timestring(3:end));
                timestringref = '1558.066';
                timediffref = str2num(timestringref(1:2))*60 + str2num(timestringref(3:end));



                timediff = timediffref - timediff; 
                lagdiff = floor(timediff*96000);

                % timeref = (0:length(ref)-1) ;
                timedata = (0:length(data)-1) - timediff;
                timestamps = timestamps - timediff;
                if timestamps(1,1) < 0; 
                    timestamps(1,1) = 1;
                end

                disp(strcat('timediff...',num2str(timediff)))
                timestamps
                % %***   DEBUG   ***%
                % % disp(strcat(num2str(size(timeref)), num2str(size(refLockout))))
                % % disp(strcat(num2str(size(timedata)), num2str(size(dataLockout))))
                
                % % timereflockout = (0:length(refLockout)-1)/fs ;
                % % timedatalockout = (0:length(dataLockout)-1)/fs + lagdiff/fs;
        
        
                % % figure(10000); grid on; hold on;
                % % plot(timedatalockout,dataLockout(:,1),'k')
                % % hold on;
                % % plot(timereflockout,refLockout(:,1),'b')
                % % title('Leadout tracks lined up')
                % % xlabel('time [s]')
                
                % %***   DEBUG ENDS  ***%
               
        
            %~~~~~~~~~~~~~~~~~~~~ NORMALIZATION ~~~~~~~~~~~~~~~~~~~~%
    
                %~~~~ separate out the 1 kHz track
                t = 1;
                fs = 96000;
                floor(timestamps(t,1)*fs)
                floor(timestamps(t,2)*fs)
                sigtime = timedata(floor(timestamps(t,1)*fs):floor(timestamps(t,2)*fs));
                sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);
        
                % sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
                % sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
                L = 2^16;
                seg = sig(floor(length(sig)/2) - L/2:floor(length(sig)/2) + L/2 - 1,:);
                
                % figure(100) 
                % plot(seg)
    
                %~~~~ window the data
                win = flattopwin(L);
                winseg = seg.*win;
                windowfactor = 0.2155774;% 0.2155774 for flattop window, 0.5 for hann window
    
                % normalization = max(abs(fftseg)); % normalize to 7 cm/s peak = 0 dB
                % normalization = max(abs(fftseg))/sqrt(2); % normalize to 5 cm/s rms = 0 dB
    
                [fftseg, fftfreq] = audio_spectrum(winseg/windowfactor, fs, 1, L);
                % figure(101)
                % audio_plotspectrum(fftfreq, fftseg, 'before normalization')
    
                w1 = 2*707/fs; w2 = 2*1404/fs;
                [b,a] = butter(4, [w1 w2]);
                segfilt = filter(b,a,seg);
                normalization = rms_response(seg);
    
    
                normalization_L = normalization(1);
                normalization_R = normalization(2);
                
                disp(strcat('max data before...', num2str(max(data))))
                disp(strcat('max data before dB...', num2str(20.0*log10(max(abs(fftseg))))))
                disp(strcat('rms data before...', num2str(rms(data))))
                disp(strcat('rms data before dB...', num2str(20.0*log10(rms(data)))))
                disp(strcat('normalization_L...', num2str(normalization_L)))
                disp(strcat('normalization_R...', num2str(normalization_R)))
                
                data(:,1)=data(:,1)./normalization_L;% now normalized to 40cm/s peak    
                data(:,2)=data(:,2)./normalization_R;% now normalized to 40cm/s peak 
    
                t = 1;
                % sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
                sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);
                L = 2^16;        
                seg = sig(floor(length(sig)/2) - L/2:floor(length(sig)/2) + L/2 - 1,:);
    
                [fftseg, fftfreq] = audio_spectrum(seg, fs, 1, L);
    
                disp(strcat('max data after...', num2str(max(data))))
                disp(strcat('max data after dB...', num2str(20.0*log10(max(abs(fftseg))))))
                disp(strcat('rms data after...', num2str(rms(data))))
                disp(strcat('rms data after dB...', num2str(20.0*log10(rms(data)))))
    
                sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);
                seg = sig(floor(length(sig)/2) - L/2:floor(length(sig)/2) + L/2 - 1,:);
                [spec, fftfreq] = audio_spectrum(seg.*win, fs, 1, 2^16);
    
                % figure(103)
                % audio_plotspectrum(fftfreq, spec, 'after normalization')
    
    
    
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
        
                    %     % floor(timestamps(t-1,1) ) - lagdiff
                    %     % floor(timestamps(t-1,2) ) - lagdiff
                    %     refT = ref(floor(timestamps(t-1,1)*fs) - lagdiff : floor(timestamps(t-1,2)*fs) - lagdiff,:);
                    % end

                    
                    if t == 1
                        sig = data(1 : floor(timestamps(1,1)*fs),:);
                        sigtime = timedata(1 : floor(timestamps(1,1)*fs)); 
                        refT = ref(1 : floor(timestamps(1,1)*fs),:);
                    elseif t == length(signal_names)
                        sig = data(floor(timestamps(end,2)*fs) : length(data),:);
                        sigtime = timedata(floor(timestamps(end,2)*fs) : length(data));  
        
                        refT = ref(floor(timestamps(end,2)*fs) : length(ref),:);
                    else
                        sig = data(floor(timestamps(t-1,1)*fs) : floor(timestamps(t-1,2)*fs),:);
                        sigtime = timedata(floor(timestamps(t-1,1)*fs) :floor(timestamps(t-1,2)*fs));  
        
                        % floor(timestamps(t-1,1) )
                        % floor(timestamps(t-1,2) )
                        refT = ref(floor(timestamps(t-1,1)*fs) : floor(timestamps(t-1,2)*fs),:);
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
                output = containers.Map(signal_names, signals)
                info_array = [lagdiff, normalization_L, normalization_R]
                disp('EXITING SEPERATE TRACKS')
    
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

