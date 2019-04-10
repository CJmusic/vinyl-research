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
    % ref_cohere = reference.tracks('transition');
    % ref_cohere = ref_cohere(coh_start*reference.fs:coh_end*reference.fs,:); 

    wave_names = [];

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
        [~,~,lagDiff] = audio_clicklineup(record.tracks('transition'), record.fs, clicks_ref);
        record.lagdiff = lagDiff;
        record.lagcorrect();

        % record.clickdetect()
        % record.clicklineup(clicks_ref)
        % record.lagcorrect()

        record.process_tracks();
        track_names = keys(record.tracks) ;
        track_data  = values(record.tracks) ;
        for j = (1:length(record.tracks));

            %% line up the audio via click detection


            %%% RMS VALUES TO CSV 
            RMS_value = rms(track_data{j});

            %%% WAVEFORM PLOTTING 
            track_time = (0:length(track_data{j})-1)/record.fs;
            figure(1); grid on; hold on;
            plot(track_time, track_data{j});
            saveas(figure(1),strcat(wave_files(i).name,'wave.png'))


            %%% SPECTRUM PLOTTING
            freq_fft = fs*(0:(n_sam/2))/n_sam;
            data_fft = fft(data(start_sam:start_sam+n_sam, :))/n_sam;
            data_fft = data_fft(1:size(data_fft)/2-1);
            plot(freq, 20.0*log10(data_fft))  
            set(gca, 'XScale', 'log');
            title(title_string)
            xlabel('Frequency (Hz)')
            ylabel('Level (dB)')  

            %%% COHERENCES TO REFERENCE RECORD
            [ amp_coh, freq_coh ] = audio_mscohere(record.tracks(tracknames{j}), reference.tracks(tracknames{j}), reference.fs);

            figure(2); grid on; hold on;
            plot(freq_coh,amp_coh(:,1))
            set(gca, 'XScale', 'log');
            xlabel('frequency [Hz]')
            title('Coherences, Left Channel')
            wave_names = [wave_names, wave_files(i).name];
            saveas(fig,filename,formattype)

            figure(3); grid on; hold on;
            plot(freq_coh,amp_coh(:,2))
            set(gca, 'XScale', 'log');
            title('Coherences, Right Channel')
            saveas(fig,filename,formattype)

            %%% COHERENCES TO PREVIOUS RECORD
            [ amp_coh, freq_coh ] = audio_mscohere(record.tracks(tracknames{j}), previous.tracks(tracknames{j}), previous.fs);

            figure(2); grid on; hold on;
            plot(freq_coh,amp_coh(:,1))
            set(gca, 'XScale', 'log');
            xlabel('frequency [Hz]')
            title('Coherences, Left Channel')
            saveas(fig,filename,formattype)

            figure(3); grid on; hold on;
            plot(freq_coh,amp_coh(:,2))
            set(gca, 'XScale', 'log');
            title('Coherences, Right Channel')
            saveas(fig,filename,formattype)
        end
        % [cdata, ctime, lagDiff] = audio_clicklineup(record.tracks('transition'), record.fs, clicks_ref); %need to make this a method of the record class

        % plot the coherence for the left and right channels compared to the reference 

        previous = record; % set the current record class to the previous one for analysis
    end

    figure(1)
    legend(wave_names)

    figure(2)
    legend(wave_names)

    figure(3)
    legend(wave_names)

end % function record_process