% This file will look at a sample recording of ANY recording and perform all the necessary measurements/analysis. 
%
% christopher zaworski
% last edit : march 31, 2019
%
%

function wav_process(folder,ref);
    clc;close all;
    disp('------------freqcohere.m---------------')
    addpath('audio_functions')
    addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');
    path_folder = strcat('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/', folder, '/')
    folder

    references = {
    '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/Bcorr/Bcorrelation_test_1.wav' 
    };    
    reference_file = references{ref}

    %addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/');

    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/');
    %record_dir = dir('');

    wave_files = dir(strcat(path_folder,'*.wav'))

    reference = audio_recordclass(reference_file)

    coh_start = 7.0;
    coh_end = 15.0;

    clicks_ref = audio_clickdetect(reference.data, reference.fs);
    disp('number of clicks in reference: ')
    size(clicks_ref)
    ref_cohere = reference.data(coh_start*reference.fs:coh_end*reference.fs,:); 
    time_ref = (0:length(ref_cohere)-1)/reference.fs; 
    % figure(1);
    % plot(time_ref, ref_cohere, 'g', 'LineWidth', 3)

    wave_names = cell(1, length(wave_files));


    for i = (1:length(wave_files)); 
        file_path = strcat(wave_files(i).folder,'/',wave_files(i).name)
        wave_names{i} = sprintf(wave_files(i).name);
        record = audio_recordclass(file_path)

        % [cdata, ctime, lagDiff] = audio_clicklineup(record.tracks('transition'), record.fs, clicks_ref); %need to make this a method of the record class
        record.clickdetect();
        disp('number of clicks: ')
        size(record.clicks)
        record.clicklineup(clicks_ref);

        % record.offset = lagDiff/record.fs + reference.offset
        % record.process_tracks();


        
        % take the proper portion of the recording to calculate the coherence 
        % the two arrays MUST be the same size
        rec_cohere = record.data;
        rec_cohere = rec_cohere(coh_start*record.fs:coh_end*record.fs,:);
        time = (0:length(rec_cohere)-1)/record.fs;

        [ amp_coh, freq_coh ] = audio_mscohere(ref_cohere, rec_cohere, reference.fs);

        % tg = uitabgroup; % opens figures in a tabbed interface 

        figure(1); hold on; grid on;
        plot(time, rec_cohere); 
        title('Records Waveforms')

        % plot the coherence for the left and right channels 
        figure(2); grid on; hold on;
        plot(freq_coh,amp_coh(:,1))
        set(gca, 'XScale', 'log');
        xlabel('frequency [Hz]')
        title('Coherences, Left Channel')

        figure(3); grid on; hold on;
        plot(freq_coh,amp_coh(:,2))
        set(gca, 'XScale', 'log');
        title('Coherences, Right Channel')

    end

    figure(1)
    legend(wave_names)

    figure(2)
    legend(wave_names)

    figure(3)
    legend(wave_names)

end % function record_process