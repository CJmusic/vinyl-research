% This file will look at a sample recording of ANY recording and perform all the necessary measurements/analysis. 
%
% christopher zaworski
% last edit : march 31, 2019
%
%

function wav_process(folder);
    clc;close all;
    disp('-----------wav_process.m---------------')
    addpath('audio_functions')
    addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');
    path_folder = strcat('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/', folder, '/')
    folder

    references = {
    '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/Bcorr/Bcorrelation_test_1.wav' 
    };    
    % reference_file = references{ref}

    %addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/');

    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/');
    %record_dir = dir('');

    wave_files = dir(strcat(path_folder,'*.wav'))


    coh_start = 7.0;
    coh_end = 15.0;
    % figure(1);
    % plot(time_ref, ref_cohere, 'g', 'LineWidth', 3)

    wave_names = cell(1, length(wave_files));


    for i = (1:length(wave_files)); 
        file_path = strcat(wave_files(i).folder,'/',wave_files(i).name)
        wave_names{i} = sprintf(wave_files(i).name);
        record = audio_recordclass(file_path)
        if i == 1; % choose first wav file as the reference to line up all the other files
            % reference = audio_recordclass(reference_file)
            reference = record;
            clicks_ref = audio_clickdetect(reference.data, reference.fs);
            disp('number of clicks in reference: ')
            size(clicks_ref)
            ref_cohere = reference.data(coh_start*reference.fs:coh_end*reference.fs,:); 
            time_ref = (0:length(ref_cohere)-1)/reference.fs; 
        end 

        %% this is to try the click lineup method
        % record.clickdetect();
        % record.clicklineup(clicks_ref);
        % disp('Click lagdiff: ')
        % figure(20); grid on; hold on;
        % plot(record.time, record.data(:,1))
        % title('Click lineup')

        % record.lagdiff = -1*record.lagdiff;
        % record.lagcorrect;

        % I need to plot what the actual correlation function looks like
        % and possibly experiment by giving it a smaller chunck of audio 
        % to calculate the correlation function with 
        xcorr_diff = audio_lineup(record.data, reference.data, record.fs);
        record.lagdiff = xcorr_diff;
        record.lagcorrect()

        % take the proper portion of the recording to calculate the coherence 
        % the two arrays MUST be the same size
        rec_cohere = record.data;
        rec_cohere = rec_cohere(coh_start*record.fs:coh_end*record.fs,:);
        time = (0:length(rec_cohere)-1)/record.fs;

        [ amp_coh, freq_coh ] = audio_mscohere(ref_cohere, rec_cohere, reference.fs);

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