%% This file runs all the calibration and testing on a recording of a test 
%  record. 
  
disp('-----------RecordProcess.m------------')
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
addpath('/Volumes/AUDIOBANK/audio_files/')
addpath('/Volumes/AUDIOBANK/audio_files/pressings/')

clc; clear all; close all; 
referencePath = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav' 
recordPath = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r30a.wav'

reference = audio_recordclass(referencePath);
reference.process_tracks()
record = audio_recordclass(recordPath);
record.process_tracks()

disp('CALLING FUNCTION')
[rmsValues, clicks] = RecordProcessTest(record, reference); 

function [rmsValues, Clicks] = RecordProcessTest(record, reference);

    disp('inside function')



    record.tracks('transition');
    reference.tracks('transition');
    class(record.tracks('leadout'))
    % data = record.tracks('leadout');
    % dataRef = reference.tracks('leadout');
    record.lagdiff = audio_corrlineup(record.tracks('leadout'), reference.tracks('leadout'), reference.fs);

    record.lagcorrect()
    record.timediff = record.lagdiff/record.fs
    record.process_tracks()

    % reference
    keys(reference.tracks)

    figure(1); hold on;
    plot(reference.track_times('1kHz'),reference.tracks('1kHz'));
    plot(record.track_times('1kHz'),record.tracks('1kHz'));
    grid on;

    figure(2); hold on;
    plot(reference.track_times('leadout'),reference.tracks('leadout'))
    plot(record.track_times('leadout'),record.tracks('leadout'))
    grid on;


    rmsValues = [];
    clicks = [];

end %RecordProcess




% function [rmsValues, Clicks] = RecordProcessTest(record)
% function [rmsValues, Clicks] = RecordProcess(record)
        %%% NORMALIZE THE AUDIO 
%         record.process_tracks;
%         data = record.tracks('1kHz');
%         level = rms(data(10*record.fs:30*record.fs)); % just take the rms level of a decent portion of the 7 cm/s sine 
    
%         trackNames = keys(record.tracks)
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
%         for i = (1:length(trackNames));
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