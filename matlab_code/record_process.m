% This file will look at a sample recording of our test record and perform all the necessary measurements/analysis. 
%
% christopher zaworski
% last edit : march 31, 2019
%
%

addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/');

addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/');
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');
record_dir = dir('');


wave_files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/A0000B0000/*a.wav'); 
wave_files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/*a.wav'); 
%this is the directory that the records are recorded under, an example is provided
% count the number of files 
%[wave_files] = dir(record_dir); 
wave_files
%reference_file = strfind(wave_files.name, '*reference.wav');
%reference_file = ['/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/reference.wav']; 

%reference_file = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/020818_A0000B0000/02072019_A0000B000r25-A.wav'; 
reference_file = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav'; 


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
coh_start = coh_start;
coh_end = coh_end;

clicks_ref = audio_clickdetect(reference.tracks('transition'), reference.fs);
ref_cohere = reference.tracks('transition');
ref_cohere = ref_cohere(coh_start*reference.fs:coh_end*reference.fs,:); 

for i = (1:length(wave_files)); 
    file_path = strcat(wave_files(i).folder,'/',wave_files(i).name)
    record = audio_recordclass(file_path)

    record.process_tracks();

    [cdata, ctime, lagDiff] = audio_clicklineup(record.tracks('transition'), record.fs, clicks_ref); %need to make this a method of the record class

    record.offset = lagDiff/record.fs + reference.offset
    record.process_tracks();

    time = (1:length(record.tracks('transition')))/record.fs;

    figure(1); hold on; grid on;
    plot(time, record.tracks('transition')); 
    title('Records 25-30, waveforms')
    
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
    title('Groove Coherences, Left Channel')

    figure(3); grid on; hold on;
    plot(freq_coh,amp_coh(:,2))
    set(gca, 'XScale', 'log');
    title('Groove Coherences, Right Channel')
end

figure(1)
legend(AUDIO_FILES)

figure(2)
legend(AUDIO_FILES)