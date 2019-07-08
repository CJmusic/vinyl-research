%% This file runs all the calibration and testing on a recording of a test 
%  record. 



function [RMS_Values, Clicks] = RecordProcess(record)

        % if CSV == true; 

        % end % if CSV true 

        % file_path = strcat(wave_files(i).folder,'/',wave_files(i).name);
        % record = audio_recordclass(file_path);

        %% line up the audio 



        %%% NORMALIZE THE AUDIO 
        data = record.tracks('1kHz');
        level = rms(record.tracks('1kHz')(10*record.fs:30*record.fs)); % just take the rms level of a decent portion of the 7 cm/s sine 
    
        record.process_tracks();
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

        for i = (1:length(record.tracks));
            track = record.tracks(i)

            rmsL = rms(track(length(track)/2 - length(track)/4:length(track)/2 + length(track)/4),1); 
            rmsR = rms(track(length(track)/2 - length(track)/4:length(track)/2 + length(track)/4),2);

            clicksL = audio_clickdetect(track(:,1)); 
            clicksR = audio_clickdetect(track(:,2)); 

            rmsValues = [rmsValues, [rmsL, rmsR]]
            Clicks = [Clicks, [clicksL, clicksR]]

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