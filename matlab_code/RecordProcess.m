%% This file runs all the calibration and testing on a recording of a test 
%  record. 
  
disp('-----------RecordProcess.m------------')
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
addpath('/Volumes/AUDIOBANK/audio_files/')
addpath('/Volumes/AUDIOBANK/audio_files/pressings/')

filepath = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r28a.wav'
record = audio_recordclass(filepath);
disp('CALLING FUNCTION')
[rmsValues, clicks] = RecordProcessTest(record); 

function [rmsValues, Clicks] = RecordProcessTest(record)
% function [rmsValues, Clicks] = RecordProcess(record)

        % if CSV == true; 

        % end % if CSV true 

        % file_path = strcat(wave_files(i).folder,'/',wave_files(i).name);
        % record = audio_recordclass(file_path);

        %% line up the audio 
        % record.process_tracks


        %%% NORMALIZE THE AUDIO 

        record.process_tracks;
        data = record.tracks('1kHz');
        level = rms(data(10*record.fs:30*record.fs)); % just take the rms level of a decent portion of the 7 cm/s sine 
    
        trackNames = keys(record.tracks)
        trackData = values(record.tracks)
        %% LOOP THROUGH EACH TRACK RECORD 

        % first the recording needs to be lined up to the reference
        % and then sliced into tracks 


        % RMS LEVEL of each track
        % Number of clicks in each track 
        % noise levels in silence transitions
        
        % COHERENCES will not be done in Record Process, they should be done 
        % seperately if needed 
        % coherences with referece and possibly subsqequent records 

        rmsValues = []; 
        Clicks = [];

        disp('BEFORE FOR LOOP')
        for i = (1:length(trackNames));
            disp('IN FOR LOOP')
            track = trackData(i);
            track = track{1};

            rmsL = rms(track(length(track)/2 - length(track)/4:length(track)/2 + length(track)/4),1); 
            rmsR = rms(track(length(track)/2 - length(track)/4:length(track)/2 + length(track)/4),2);

            %% A WEIGHTING AND CCIR WEIGHTING RMS VALUES 

            clicksL = audio_clickdetect(track(:,1), record.fs); 
            clicksR = audio_clickdetect(track(:,2), record.fs); 

            %% DIVIDE SPECTRUM INTO CHUNKS (ARM RESONANCE, HIGH COHERENCE (1000-300Hz) and HIGH FREQ)
            %  SAVE dB power in these regions

            rmsValues = [rmsValues, [rmsL, rmsR]];
            Clicks = [Clicks, [clicksL, clicksR]];


            % if CSV = true;
            %     % WRITE ALL THE APPROPRIATE VALUES TO THE CSV FILE 
            % end % if CSV = true

            % coarseness = 20; % number of samples we allow the mode to be off by when counting common clicks 
            % [dt_row, dt_column] = find(click_matrix > (lagdiff - coarseness) & click_matrix < (lagdiff + coarseness));
            % clicks_total = size(click_matrix,2);
            % clicks_common = length(dt_row);
            % clicks_unique = length(dt_row) -  size(clicks, 2);  

            % RMS_valueL = rms(record.tracks('transition'));
            % RMS_valueR = rms(record.tracks('transition'));

            % CSV_MATRIX = [CSV_MATRIX ; {record.filename, 'transition', RMS_valueL, RMS_valueR, clicks_total, clicks_common, clicks_unique}];

        end % for track loop 

        %% IF Individual TRACK 
        %%% FOR EACH INDIVIDUAL TRACK DESIGN A SERIES OF 
        %   NOTCH FILTERS to measure background noise
        %   this will tell us the signal to noise of the 
        %   test disk 


        %% SIGNAL TO NOISE FOR EACH TRACKS
        % measure the strength of the harmonics by counting bins that contain signal
        % and bins that contain noise, dividing them by number of bins and comparing

        harmonics = (1:20)/2 % 20 harmonics ?  
        n_sam = 2^16;
        fs = record.fs;

        %% 1 kHz
        disp('at 1kHz')
        data = record.tracks('1kHz');
        data_fft = fft(data)/n_sam;
        data_fft = data_fft(1:n_sam/2+1);
        % 1k_bins =(harmonics*1000/fs); % use the floor operator to get integer bins
        kbins = floor(harmonics*1000)
        n_sam = 2^16;
        freq_fft = fs*(0:(n_sam/2))/n_sam;

        k = sum(data_fft(kbins))/length(harmonics);
        data_fft(kbins) = 0.0;
        noise = sum(data_fft)/(length(data_fft) - length(harmonics));

        % for i = (1:length(1k))



        %% 10 kHz 


        %% 100 Hz


        %% 3150 Hz 






        % CSV_TABLE = cell2table(CSV_MATRIX)
        % writetable(CSV_TABLE,strcat(reference.directory,'A0000B0000_analysis.csv'))
        % figure(30)
        % legend(wave_files.name)
        % figure(40)
        % legend(wave_files.name)
        % saveas(figure(30),strcat(record.directory, '/plots/','Cumulative_cohL.png'))
        % saveas(figure(40),strcat(record.directory, '/plots/','Cumulative_cohR.png'))


        %%IF LOG SWEEP 
        % freqResp = fft(record.track(SWEEP))

end %record process