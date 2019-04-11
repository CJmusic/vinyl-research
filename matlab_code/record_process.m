% This file will look at a sample recording of our test record and perform all the necessary measurements/analysis. 
%
% christopher zaworski
% last edit : march 31, 2019
%
%

function record_process(folder,ref);
    clc;close all;
    disp('-----------recordprocess.m---------------')
    addpath('audio_functions')
    addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');
    path_folder = strcat('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/', folder)
    folder
    references = {'/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav'};    
    reference_file = references{ref}

    %addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/');

    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/');
    %record_dir = dir('');

    if strcmp(folder,'A0000B0000');
        wave_files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/*a.wav'); 
        % wave_files = dir('/Volumes/AUDIOBANK/audio_files/A0000B0000/*a.wav');
    else
        wave_files = dir(strcat(path_folder,'*.wav')); 
    end
    wave_files

    % string the file names of relevant data 

    %process the string and pull info: 
    pressing_date = 0; 
    recording_date = 0; 
    recording_timestamp = 0; 

    top_stamper = 0;
    top_stamper_hits = 0;
    bottom_stamper = 0;
    bottom_stamper_hits = 0;

    %locate in the csv file

    path_csv = 0; 
    csv_file = 0; 

    reference = audio_refrecordclass(reference_file)
    %signal_array = audio_detectsignal(reference.dataL); 

    % figure out where the signals are in the reference track
    % line up all signals with the reference
    % use the reference ref_timestamps

    % I also need to add the second set of signals (since they repeat, and the ref_transition)
    % for some reason theres a second log sweep

    % take the clicks array and the portion of the record used to calculate the coherence
    % coh_start = coh_start;    
    % coh_end = coh_end;

    clicks_ref = audio_clickdetect(reference.tracks('transition'), reference.fs);

    wave_names = [];

    CSV_titles = {'wavefile','track name', 'RMS', 'total clicks', 'common clicks', 'unique clicks'};
    CSV_MATRIX = cell(1,6);
    CSV_MATRIX(1,:) = CSV_titles;

    for i = (1:length(wave_files)); 
        disp('~~~~~~~~~~NEXT FILE~~~~~~~~~~~~')
        file_path = strcat(wave_files(i).folder,'/',wave_files(i).name)
        wave_names = [wave_names, wave_files(i).name];
        disp('wave_name')
        wave_name = wave_files(i).name(1:end-4) 
        record = audio_recordclass(file_path)
        if i == 1; % choose first wav file as the reference to line up all the other files
            previous = record;
        end 

        record.transition_track()
        disp('Transition track')
        size(record.tracks('transition'))
        clicks_tran = audio_clickdetect(record.tracks('transition'), record.fs);
        % [lagDiff] = audio_clicklineup(record.tracks('transition'), record.fs, clicks_ref);
        [click_matrix_tran, lagDiff] = audio_clickmatrix(clicks_tran, clicks_ref);


        record.lagdiff = lagDiff;
        record.lagcorrect()
        record.process_tracks();

        track_names = keys(record.tracks) ;
        track_data  = values(record.tracks) ;

        for j = (1:length(record.tracks));
            %%% RMS VALUES TO CSV 
            RMS_value = rms(track_data{j});
            % RMS_values = [RMS_values, RMS_value];

            %%% WAVEFORM PLOTTING 
            track_time = (0:length(track_data{j})-1)/record.fs;
            figure(1); grid on; 
            plot(track_time, track_data{j});
            title(strcat(wave_name,track_names{j},'Waveform'))
            xlabel('time (s)')
            ylabel('Amplitude')  
            saveas(figure(1),strcat(wave_name,track_names{j},'wave.png'))


            %%% SPECTRUM PLOTTING
            disp('SPECTRUM PLOTTING')
            n_sam = length(track_data{j})
            
            freq_fft = record.fs*(0:(n_sam/2))/n_sam;
            
            data_fft = fft(track_data{j})/n_sam;
            size(freq_fft)
 
            data_fft = data_fft(1:n_sam/2+1);
            size(data_fft)
            figure(2); grid on; 
            plot(freq_fft, 20.0*log10(data_fft))  
            set(gca, 'XScale', 'log');
            title(strcat(wave_name,track_names{j},'Spectrum'))
            xlabel('Frequency (Hz)')
            ylabel('Level (dB)')  
            saveas(figure(2),strcat(wave_name,track_names{j},'spectrum.png'))

            %%% COHERENCES TO REFERENCE RECORD
            [ amp_coh, freq_coh ] = audio_mscohere(record.tracks(track_names{j}), reference.tracks(track_names{j}), reference.fs);

            figure(3); grid on; 
            plot(freq_coh,amp_coh(:,1))
            set(gca, 'XScale', 'log');
            xlabel('frequency [Hz]')
            % title('Coherences to Reference, Left Channel')
            title(strcat(wave_name,track_names{j},'Coherences to Reference, Left Channel'))
            wave_names = [wave_names, wave_files(i).name];
            filename = strcat(wave_name,track_names{j},'cohleft_ref.png')
            saveas(figure(3),filename)

            figure(4); grid on; 
            plot(freq_coh,amp_coh(:,2))
            set(gca, 'XScale', 'log');
            % title('Coherences to Reference, Right Channel')
            title(strcat(wave_name,track_names{j},'Coherences to Reference, Right Channel'))
            filename = strcat(wave_name,track_names{j},'cohright_ref.png')
            saveas(figure(4),filename)

            %%% COHERENCES TO PREVIOUS RECORD
            [ amp_coh, freq_coh ] = audio_mscohere(record.tracks(track_names{j}), previous.tracks(track_names{j}), previous.fs);

            figure(5); grid on; 
            plot(freq_coh,amp_coh(:,1))
            set(gca, 'XScale', 'log');
            xlabel('frequency [Hz]')
            % title('Coherences to Previous, Left Channel')
            title(strcat(wave_name,track_names{j},'Coherences to Previous, Right Channel'))
            filename = strcat(wave_name,track_names{j},'cohleft_prev.png')
            saveas(figure(5),filename)

            figure(6); grid on; 
            plot(freq_coh,amp_coh(:,2))
            set(gca, 'XScale', 'log');
            title('Coherences to Previous, Right Channel')
            title(strcat(wave_name,track_names{j},'Coherences to Previous, Right Channel'))
            filename = strcat(wave_name,track_names{j},'cohleft_prev.png')
            saveas(figure(6),filename)

            if strcmp(track_names{j},'transition');
                clicks = audio_clickdetect(record.tracks('transition'), record.fs);
                [click_matrix, lagdiff] = audio_clickmatrix();

                coarseness = 20; % number of samples we allow the mode to be off by when counting common clicks 
                [dt_row, dt_column] = find(click_matrix > lagDiff - coarseness && click_matrix < lagDiff + coarseness);
                clicks_total = size(click_matrix,1);
                clicks_common = len(dt_row);
                clicks_unique = size(clicks, 1) - len(dt_row);

                CSV_MATRIX = {CSV_MATRIX ; strcat(wave_name), track_names{j}, RMS_value, clicks_total, clicks_common, clicks_unique};

            else; 
                track_name = track_names{j}
                CSV_MATRIX(end+1,:) = {wave_name, track_names, RMS_value, 'n/a', 'n/a', 'n/a'};
            
            end 
            clf(figure(1))
            clf(figure(2))
            clf(figure(3))
            clf(figure(4))
            clf(figure(5))
            clf(figure(6))
        end % for loop through record tracks

        previous = record; % set the current record class to the previous one for analysis
    end % for loop of record files

    % Write the table to a CSV file
    CSV_TABLE = cell2table(CSV_MATRIX,'VariableNames',CSV_titles)
    writetable(CSV_TABLE,'A0000B0000_analysis.csv')

end % function record_process