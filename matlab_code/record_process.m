% This file will look at a sample recording of our test record and perform all the necessary measurements/analysis. 
% Christopher Zaworski
% Last edit : March 7, 2019
%
%

addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/');

record_dir = dir('');


wave_files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/*.wav'); 
%this is the directory that the records are recorded under, an example is provided
% count the number of files 
%[wave_files] = dir(record_dir); 
wave_files
%reference_file = strfind(wave_files.name, '*reference.wav');
%reference_file = ['/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/reference.wav']; 

reference_file = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/020818_A0000B0000/02072019_A0000B000r25-A.wav'; 


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


% looks like ref_signals is 1x26 
% and signals is 1x15



%ref_signals  = containers.map(signals(i), [ref_signalsL, ref_signalsR]; % probably wont work, I need a list of these
                                                                        % arrays above and map keys to them all at once

%%% This code left for reference, it should plot out each segement that's been cut from the timestamps and track listing
%for i = (1:length(ref_timestamps)-1);
%    data_seg = reference.dataL(floor(ref_timestamps(i)*reference.fs):floor(ref_timestamps(i+1)*reference.fs));
%    time_seg = time(floor(ref_timestamps(i)*reference.fs):floor(ref_timestamps(i+1)*reference.fs));
%    plot(time_seg, data_seg);
%    line([ref_timestamps(i) ref_timestamps(i)], get(gca, 'ylim'),'Color', 'blue','LineStyle', '--');
%    
%    data_seg = reference.dataL(floor(ref_timestamps2(i)*reference.fs):floor(ref_timestamps2(i+1)*reference.fs));
%    time_seg = time(floor(ref_timestamps2(i)*reference.fs):floor(ref_timestamps2(i+1)*reference.fs));
%    plot(time_seg, data_seg);
%    line([ref_timestamps2(i) ref_timestamps2(i)], get(gca, 'ylim'),'Color', 'blue','LineStyle', '--');
%end


% figure out where the signals are in the reference track
% line up all signals with the reference
% use the reference ref_timestamps
%
% I also need to add the second set of signals (since they repeat, and the ref_transition)
% for some reason theres a second log sweep

%for i = (1:length(wave_files)); 
%    file_path = strcat(wave_files(i).folder,'/',wave_files(i).name)
%    record = audio_recordclass(file_path)
%    signal_array = audio_detectsignal(record.data); 
%    %click_mat = audio_clickmatrix(record.clicks, reference.clicks)
%    
%    %record_id = 'r*.wav';

%    %[data, fs] = audioread(file_path); 
%    %[clicks] = audio_clickdetect(data, fs);
%    %lagdiff = mode(click_mat);
%    %click_info = audio_clickcompare(clicks, clicks) 
%    %[cdata, ctime, cdata_ref, ctime_ref] = audio_lineup(data, fs, data_ref, fs_ref, lagdiff);
%end


