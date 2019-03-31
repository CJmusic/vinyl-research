% This file will look at a sample recording of any recording and perform all the necessary measurements/analysis. 
%
% christopher zaworski
% last edit : march 31, 2019
%
%

function wav_process(folder,ref);
    addpath('audio_functions')
    addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');
    path_folder = strcat('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/', folder)
    folder

    references = {'/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav',
    '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/Bcorr/Bcorrelation_test_1.wav' 
    };    
    reference_file = references{ref}

    %addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/');

    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/');
    %record_dir = dir('');

    wave_files = dir(strcat(path_folder,'*.wav'))

    reference = audio_refrecordclass(reference_file)

    coh_start = coh_start;
    coh_end = coh_end;

    clicks_ref = audio_clickdetect(reference.tracks('transition'), reference.fs);
    ref_cohere = reference.tracks('transition');
    ref_cohere = ref_cohere(coh_start*reference.fs:coh_end*reference.fs,:); 

    wave_names = [];

    for i = (1:length(wave_files)); 
        file_path = strcat(wave_files(i).folder,'/',wave_files(i).name)
        wave_names = [wave_names, wave_files(i).name];
        record = audio_recordclass(file_path)

        record.process_tracks();

        [cdata, ctime, lagDiff] = audio_clicklineup(record.tracks('transition'), record.fs, clicks_ref); %need to make this a method of the record class

        record.offset = lagDiff/record.fs + reference.offset
        record.process_tracks();

        time = (1:length(record.tracks('transition')))/record.fs;

        figure(1); hold on; grid on;
        plot(time, record.tracks('transition')); 
        title('Records Waveforms')
        
        % take the proper portion of the recording to calculate the coherence 
        % the two arrays MUST be the same size
        rec_cohere = record.tracks('transition');
        rec_cohere = rec_cohere(coh_start*record.fs:coh_end*record.fs,:);
        [ amp_coh, freq_coh ] = audio_mscohere(ref_cohere, rec_cohere, reference.fs);

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