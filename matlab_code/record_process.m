% This file will look at a sample recording of our test record and perform all the necessary measurements/analysis. 
% Christopher Zaworski
% Last edit : March 7, 2019
%
%

addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/');

record_dir = dir('');


wave_files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/A0000B0000/*a.wav'); 
%this is the directory that the records are recorded under, an example is provided
% count the number of files 
%[wave_files] = dir(record_dir); 
wave_files
%reference_file = strfind(wave_files.name, '*reference.wav');
%reference_file = ['/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/reference.wav']; 

%reference_file = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/020818_A0000B0000/02072019_A0000B000r25-A.wav'; 
reference_file = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/A0000B0000/031418_A0000B0000r27a.wav'; 


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
%
% I also need to add the second set of signals (since they repeat, and the ref_transition)
% for some reason theres a second log sweep
clicks_ref = audio_clickdetect(reference.tracks('transition'), reference.fs);
for i = (1:length(wave_files)); 
    file_path = strcat(wave_files(i).folder,'/',wave_files(i).name)
    record = audio_recordclass(file_path)
    %signal_array = audio_detectsignal(record.data); 
    %click_mat = audio_clickmatrix(record.clicks, reference.clicks)
    
    %lagDiffL = audio_lineup(record.dataL(2*record.fs:5*record.fs), record.fs, reference.dataL(2*record.fs:5*record.fs))
    %lagDiffR = audio_lineup(record.dataR(2*record.fs:5*record.fs), record.fs, reference.dataR(2*record.fs:5*record.fs))
    %record.lagcorrect(lagDiffL, lagDiffR);
    
    %[cdataL, ctimeL, cdata_refL, ctime_refL] = audio_lagcorrect(record.dataL, record.fs, reference.dataL, reference.fs, lagDiffL);
    %[cdataR, ctimeR, cdata_refR, ctime_refR] = audio_lagcorrect(record.dataR, record.fs, reference.dataR, reference.fs, lagDiffR);
    disp('FIGURE TIME')
    figure(1); hold on; grid on;
   % plot(record.dataL(1:10:end))
    %plot(reference.dataL,'g')



    %record = record.process_tracks();
    record.process_tracks();
   % figure(i);hold on;
    disp('in script printing tracks')
    record.tracks
    %keys(record.tracks)
    %plot(record.tracks('transition'));
    %plot(reference.tracks('1kHz'));

    [cdata, ctime] = audio_clicklineup(record.tracks('transition'), record.fs, clicks_ref);
    plot(ctime,cdata) 
    %lagdiffL = audio_clicklineup(record.clicksL(1:100), reference.clicksL(1:100))
    %lagdiffR = audio_clicklineup(record.clicksR(1:100), reference.clicksR(1:100))
    %record_id = 'r*.wav';

    %[data, fs] = audioread(file_path); 
    %[clicks] = audio_clickdetect(data, fs);
    %lagdiff = mode(click_mat);
    %click_info = audio_clickcompare(clicks, clicks) 
    %[cdata, ctime, cdata_ref, ctime_ref] = audio_lineup(data, fs, data_ref, fs_ref, lagdiff);
end


