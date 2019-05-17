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
    % path_folder = strcat('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/', folder, '/')
    path_folder = strcat('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/', folder, '/')

    % references = {
    % '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/Bcorr/Bcorrelation_test_1.wav',
    % '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/030419_r26fivetrials/one.wav' 
    % };    
    % reference_file = references{ref}

    addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/');

    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/');
    %record_dir = dir('');

    wave_files = dir(strcat(path_folder,'*.wav'));


    coh_start = 8.0;
    coh_end = 18.0;
    % figure(1);
    % plot(time_ref, ref_cohere, 'g', 'LineWidth', 3)

    wave_names = cell(1, length(wave_files))
    for i = (1:1)%length(wave_files))
        file_path = strcat(wave_files(i).folder,'/',wave_files(i).name)
        wave_names{i} = sprintf(wave_files(i).name)
        record = audio_recordclass(file_path)
        if i == 1; % choose first wav file as the reference to line up all the other files
            reference = record;
            % clicks_ref = audio_clickdetect(reference.data, reference.fs);
        %     clicks_ref = aud_clickremoval(reference.data, reference.fs);
            ref_cohere = reference.data(coh_start*reference.fs:coh_end*reference.fs,:); 
            time_ref = (0:length(ref_cohere)-1)/reference.fs; 
        end 

        disp('BUFFER INFO')
        [data_aud, clicks] = aud_clickremoval(record.data(:,1));
        size(data_aud)
        
        figure(1)
        plot(record.time, record.data(:,1))
        grid on;

        figure(2)
        plot(record.time, data_aud)
        grid on;
        % SP=1; %your point goes here 
        % hax=axes; 
        % for j = (1:length(clicks))
        %     line([clicks(j) clicks(j)],get(hax,'YLim'),'Color',[1 0 0])
        % end


%~~~~~~~~~~~~~~~~~NORMALIZATION TEST~~~~~~~~~~~~~~~~~~~~~

        % amplitude = audio_findamplitude(record.data,1000,record.fs)
        % amplitude

%~~~~~~~~~~~~~~~NORMALIZATION TEST END~~~~~~~~~~~~~~~~~~~




%
% %~~~~~~~~~~~~~~~~~~CLICK LINEUP TEST~~~~~~~~~~~~~~~~~~~~~

%         % figure(20); grid on; hold on;
%         % plot(record.time, record.data(:,1))
%         % title('Click lineup')

%         record.lagdiff = -1*record.lagdiff;
%         record.lagcorrect;
        
%         figure(10); hold on; grid on; 
%         title('pre lineup')
%         plot(record.time, record.data)


%         % record.clickdetect();
%         % record.clicklineup(clicks_ref);
%         % record.lagcorrect()
%         % record.time = record.time - record.lagdiff/record.fs;


%         % figure(20); hold on; grid on;
%         % title('post lineup clicks')
%         % plot(record.time,record.data)

%         % record.time = record.time + record.lagdiff/record.fs;
%         % record.lagdiff = -1.0*record.lagdiff;
%         % record.lagcorrect()

%         xcorr_diff = audio_lineup(record.data, reference.data, record.fs);
%         record.lagdiff = xcorr_diff;
%         record.lagcorrect()
%         % record.time = record.time - record.lagdiff/record.fs;

%         figure(30); hold on; grid on;
%         title('post lineup xcorr')
%         plot(record.time,record.data)

% ~~~~~~~~~~~~~~~CLICK LINEUP TEST END~~~~~~~~~~~~~~~~~~~
%}


%~~~~~~~~~~~~~~~~~~~~COHERENCE TEST~~~~~~~~~~~~~~~~~~~~~~
        % take the proper portion of the recording to calculate the coherence 
%         rec_cohere = record.data;
%         rec_cohere = rec_cohere(coh_start*record.fs:coh_end*record.fs,:);
%         % time = (0:length(rec_cohere)-1)/record.fs;

%         [ amp_coh, freq_coh ] = audio_mscohere(ref_cohere, rec_cohere, reference.fs);

%         % figure(1); hold on; grid on;
%         % plot(time, rec_cohere); 
%         % title('Records Waveforms')
%         n_sam = length(rec_cohere)

%         freq_fft = record.fs*(0:(n_sam/2))/n_sam;

%         data_fft = fft(rec_cohere)/n_sam;
%         size(freq_fft)

%         data_fft = data_fft(1:n_sam/2+1);
%         size(data_fft)

%         figure(4); grid on; hold on;
%         plot(freq_fft, 20.0*log10(data_fft))  
%         set(gca, 'XScale', 'log');
%         title('Spectrum of noise')
%         xlabel('Frequency (Hz)')
%         ylabel('Level (dB)')  
%         % % plot the coherence for the left and right channels 
%         figure(2); grid on; hold on;
%         plot(freq_coh,amp_coh(:,1))
%         set(gca, 'XScale', 'log');
%         xlabel('frequency [Hz]')
%         title('Coherences, Left Channel')

%         figure(3); grid on; hold on;
%         plot(freq_coh,amp_coh(:,2))
%         set(gca, 'XScale', 'log');
%         title('Coherences, Right Channel')

% % ~~~~~~~~~~~~~~~COHERENCE TEST END~~~~~~~~~~~~~~~~~~~
%}
    end % wav files loop

%     figure(10)
%     legend(wave_names)

%     figure(20)
%     legend(wave_names)

%     figure(30)
%     legend(wave_names)

%     % figure(1)
%     % legend(wave_names)

%     figure(2)
%     legend(wave_names)

%     figure(3)
%     legend(wave_names)

end % function record_process