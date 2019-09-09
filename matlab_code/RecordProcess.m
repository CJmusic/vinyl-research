%% This file runs all the calibration and testing on a recording of a test 
%  record. 
  
disp('-----------RecordProcess.m------------')
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
addpath('/Volumes/AUDIOBANK/audio_files/')
addpath('/Volumes/AUDIOBANK/audio_files/pressings/')

clc; close all;
referencePath = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r28a.wav' 
recordPath = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r30a.wav'

reference = audio_recordclass(referencePath);
reference.process_tracks()
record = audio_recordclass(recordPath);
record.process_tracks()

disp('CALLING FUNCTION')
[rmsValues, clicks] = RecordProcessTest(record, reference); 

function [rmsValues, Clicks] = RecordProcessTest(record, reference);

        disp('inside function')

        Clicks = [];
        rmsValues = [];

        %%%~~~~~~~ LINEUP AUDIO ~~~~~~~%%%

        %% ISOLATING LINEUP %% 
        data = record.tracks('leadout');
        dataRef = reference.tracks('leadout');

        lagdiff = audio_corrlineup(data, dataRef);
        data2corr = circshift(data, lagdiff);

        time = (0:(length(data2corr)-1))/record.fs;
        timeRef = (0:(length(dataRef)-1))/record.fs;
        figure(30)
        hold on; grid on; 
        plot(timeRef, dataRef)
        plot(time, data2corr)
        %%%   ISOLATING END    %%%

        figure(1); hold on;
        title('Pre-lineup')
        plot(reference.track_times('leadout'),reference.tracks('leadout'))
        plot(record.track_times('leadout'),record.tracks('leadout'))
        grid on;

        record.lagdiff = audio_corrlineup(record.tracks('leadout'), reference.tracks('leadout'));

        % record.lagcorrect()

        % record.tracks('leadout') = circshift(record.tracks('leadout'),record.lagdiff);
        data = record.tracks('leadout');
        dataRef = reference.tracks('leadout');


        figure(2); hold on;
        title('Post-lineup')
        plot(reference.track_times('leadout'),reference.tracks('leadout'))
        plot(record.track_times('leadout'),record.tracks('leadout'))
        grid on;

        figure(3); hold on;
        plot(reference.track_times('1kHz'),reference.tracks('1kHz'));
        plot(record.track_times('1kHz'),record.tracks('1kHz'));
        grid on;


        % figure(3); hold on; 
        % % plot(data)
        % plot(dataRef)
        % plot(data2corr)
        % grid on;

        %%%~~~~~~~ LINEUP END ~~~~~~~%%%

        %%%~~~~~~~ NORMALIZE ~~~~~~~~%%%
        data = record.tracks('1kHz');
        tstart = 15;
        tend = 30;

        dataL = data(tstart*record.fs:tend*record.fs,1);
        dataR = data(tstart*record.fs:tend*record.fs,2);

        RMSdataR = sqrt(sum(dataR.^2)/((tstart-tend)*record.fs))
        RMSdataL = sqrt(sum(dataL.^2)/((tstart-tend)*record.fs))

        PEAKdataR=sqrt(2)*RMSdataR*40/7; %digital value of peak level
        PEAKdataL=sqrt(2)*RMSdataL*40/7; %digital value of peak level
        record.data = [record.data(:,1)/PEAKdataR, record.data(:,2)/PEAKdataL];

        Normalization = [PEAKdataL, PEAKdataR];

        %%%~~~~~ NORMALIZE END ~~~~~~%%%

        %%%~~~~~~ CLICK DETECT ~~~~~~%%%
        % [ClicksL, ClicksR] = length(audio_clickdetect(record.data))
        %%%~~~~ CLICK DETECT END ~~~~~~%%%
        
        %%~~~ TRACK SPECIFIC TESTS ~~~%%%
        % signal_names = keys(record.tracks);

        for i = (1:length(record.signal_names));
                % if record.signal_names(i) == ('1kHz') || record.signal_names(i) == ('1kHz2') || record.signal_names(i) == ('1kHzV') 
                if ismember(record.signal_names(i),{'1kHz','1kHzL', '1kHzR', '1kHzV','1kHz2','1kHzL2', '1kHzR2', '1kHzV2'});
                        disp('1kHz')
                        %% Flattop window to measure the peaks 
                end % 1 kHz tests


                if ismember(record.signal_names(i),{'10kHz','10kHz', '10kHzL', '10kHzR', '10kHz2','10kHz2', '10kHzL2', '10kHzR2'});
                        disp('10kHz')

                end % 10 kHz tests


                if ismember(record.signal_names(i),{'100Hz', '100HzL', '100HzR','100Hz2', '100HzL2', '100HzR2'});
                        disp('100Hz')

                end % 100 Hz tests


                if ismember(record.signal_names(i),{'sweep', 'sweepL', 'sweepR', 'sweepV', 'sweep2', 'sweepL2', 'sweepR2', 'sweepV2'});
                        disp('sweep')

                end % sweep tests
                

                if ismember(record.signal_names(i),{'3150Hz','3150Hz2'});
                        disp('3150Hz')
                end % 3150Hz tests


                if ismember(record.signal_names(i),{'transition'});
                        disp('transition')

                end % transition tests

        end % for loop signal_names 

        %%~~~ TRACK SPECIFIC TESTS ~~~%%%
end %RecordProcess




% function [rmsValues, Clicks] = RecordProcessTest(record)
% function [rmsValues, Clicks] = RecordProcess(record)
        %%% NORMALIZE THE AUDIO 
%         record.process_tracks;
%         data = record.tracks('1kHz');
%         level = rms(data(10*record.fs:30*record.fs)); % just take the rms level of a decent portion of the 7 cm/s sine 
    
%         signal_names = keys(record.tracks)
%         trackData = values(record.tracks)
%         %% LOOP THROUGH EACH TRACK RECORD 

%         % first the recording needs to be lined up to the reference
%         % and then sliced into tracks 


%         % RMS LEVEL of each track
%         % Number of clicks in each track 
%         % noise levels in silence transitions
        
%         % COHERENCES will not be done in Record Process, they should be done 
%         % seperately if needed 
%         % coherences with referece and possibly subsqequent records 

%         rmsValues = []; 
%         Clicks = [];

%         disp('BEFORE FOR LOOP')
%         for i = (1:length(signal_names));
%             disp('IN FOR LOOP')
%             track = trackData(i);
%             track = track{1};

%             rmsL = rms(track(length(track)/2 - length(track)/4:length(track)/2 + length(track)/4),1); 
%             rmsR = rms(track(length(track)/2 - length(track)/4:length(track)/2 + length(track)/4),2);

%             %% A WEIGHTING AND CCIR WEIGHTING RMS VALUES 

%             clicksL = audio_clickdetect(track(:,1), record.fs); 
%             clicksR = audio_clickdetect(track(:,2), record.fs); 

%             %% DIVIDE SPECTRUM INTO CHUNKS (ARM RESONANCE, HIGH COHERENCE (1000-300Hz) and HIGH FREQ)
%             %  SAVE dB power in these regions

%             rmsValues = [rmsValues, [rmsL, rmsR]];
%             Clicks = [Clicks, [clicksL, clicksR]];


%             % if CSV = true;
%             %     % WRITE ALL THE APPROPRIATE VALUES TO THE CSV FILE 
%             % end % if CSV = true

%             % coarseness = 20; % number of samples we allow the mode to be off by when counting common clicks 
%             % [dt_row, dt_column] = find(click_matrix > (lagdiff - coarseness) & click_matrix < (lagdiff + coarseness));
%             % clicks_total = size(click_matrix,2);
%             % clicks_common = length(dt_row);
%             % clicks_unique = length(dt_row) -  size(clicks, 2);  

%             % RMS_valueL = rms(record.tracks('transition'));
%             % RMS_valueR = rms(record.tracks('transition'));

%             % CSV_MATRIX = [CSV_MATRIX ; {record.filename, 'transition', RMS_valueL, RMS_valueR, clicks_total, clicks_common, clicks_unique}];

%         end % for track loop 

%         %% IF Individual TRACK 
%         %%% FOR EACH INDIVIDUAL TRACK DESIGN A SERIES OF 
%         %   NOTCH FILTERS to measure background noise
%         %   this will tell us the signal to noise of the 
%         %   test disk 


%         harmonics = (1:20)/2 % 20 harmonics ?  
%         n_sam = 2^16;
%         fs = record.fs;
%         data = record.tracks('1kHz');
%         data_fft = fft(data)/n_sam;
%         data_fft = data_fft(1:n_sam/2+1);
%         n_sam = 2^16;
%         freq_fft = fs*(0:(n_sam/2))/n_sam;
%         %% SIGNAL TO NOISE FOR EACH TRACKS
%         % measure the strength of the harmonics by counting bins that contain signal
%         % and bins that contain noise, dividing them by number of bins and comparing
%         signalBins = [681,682,683,684,685]

%         noiseBins = [1:length(data_fft)];

%         % noiseBins(signalBins) = [];
        
%         % disp('NoiseBins')
%         % noiseBins(signalBins);
%         % noiseBins(signalBins) = [];
%         signal = sum(data(signalBins))/length(signalBins);
%         noise = sum(data(noiseBins))/length(noiseBins);


%         %%% DISTORTION/HARMONICS measurement 
%         % 1k_bins =(harmonics*1000/fs); % use the floor operator to get integer bins
%         kbins = floor(harmonics*1000)
%         %% 1 kHz
%         disp('at 1kHz')

%         k = sum(data_fft(kbins))/length(harmonics);
%         % data_fft(kbins) = 0.0;
%         % noise = sum(data_fft)/(length(data_fft) - length(harmonics));



%         figure(1); grid on;
%         plot(data)

%         figure(10); grid on; hold on; 
%         plot(freq_fft,20*log10(real(data_fft)),'blue');
%         grid on; hold on; 
%         set(gca, 'XScale', 'log');
%         disp('plotting')
%         for xi = (1:length(signalBins))
%             x1 = freq_fft(signalBins(xi))
%             line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
%         end

%         for xi = (1:length(kbins))
%             x1 = freq_fft(kbins(xi))
%             line([x1 x1], get(gca, 'ylim'),'Color', 'green','LineStyle', '--');
%         end

%         %% 10 kHz 


%         %% 100 Hz


%         %% 3150 Hz 






%         % CSV_TABLE = cell2table(CSV_MATRIX)
%         % writetable(CSV_TABLE,strcat(reference.directory,'A0000B0000_analysis.csv'))
%         % figure(30)
%         % legend(wave_files.name)
%         % figure(40)
%         % legend(wave_files.name)
%         % saveas(figure(30),strcat(record.directory, '/plots/','Cumulative_cohL.png'))
%         % saveas(figure(40),strcat(record.directory, '/plots/','Cumulative_cohR.png'))


%         %%IF LOG SWEEP 
%         % freqResp = fft(record.track(SWEEP))

% end %record process