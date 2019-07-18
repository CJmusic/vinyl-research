%% This file processes a large number of records from a test pressing, 
%   
% 



function BatchProcess(folder)
    % clc;close all;
    disp('-----------BatchProcess.m------------')
    addpath('audio_functions')
    addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
    addpath('/Volumes/AUDIOBANK/audio_files/')
    addpath('/Volumes/AUDIOBANK/audio_files/pressings/')

    path_folder = strcat('/Volumes/AUDIOBANK/audio_files/pressings/', folder, '/')

    wave_files = dir(strcat(path_folder,'/*.wav'));

    reference = audio_recordclass('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r28a.wav');

    reference.process_tracks();

    % clicks_ref = audio_clickdetect(reference.tracks('transition'), reference.fs);

    % data_ref = reference.tracks('transition');
    % data_refL = data_ref(:,1);
    % data_refR = data_ref(:,2);

    % CSV_titles = {"wavefile", "track_name", "RMS_L", "RMS_R", "total_clicks", "common_clicks", "unique_clicks"};
    CSV_titles = {"Filename", "Pressing", "Record", "Side", "Normalization", "Track",  "RMSL", "RMSR", "ClicksL", "ClicksR"};
    CSV_MATRIX = cell(1,10);
    CSV_MATRIX(1,:) = CSV_titles;


    for i = (1:length(wave_files))
        disp('FILE PATH ')
        filepath = strcat(wave_files(i).folder,'/',wave_files(i).name)
        wave_files(i).folder
        filename = wave_files(i).name
        filepath =strcat(wave_files(i).folder,'/', filename(3:end)) 

        record = audio_recordclass(filepath);
        % record = audio_recordclass(strcat(wave_files(i).folder,'/',wave_files(i).name));

        %% GO BACK TO CLICK LINE UP 
        record.lagdiff = audio_corrlineup(record.data(1:60*record.fs), reference.tracks('1kHz'), record.fs);
        record.lagcorrect();
        % lagdiff(reference, record)
        [RMS_values, clicks] = RecordProcess(record, reference); 
    end


    CSV_TABLE = cell2table(CSV_MATRIX)
    writetable(CSV_TABLE,strcat(reference.directory,'A0000B0000_analysis.csv'))

end % function BatchProcess