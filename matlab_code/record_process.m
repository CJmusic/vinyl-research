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

reference_file = ['/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/020818_A0000B0000/02072019_A0000B000r25-A.wav']; 


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

reference = audio_recordclass(reference_file)
signal_array = audio_detectsignal(reference.dataL); 

%%% RECORD INFO 
%signal_names = {'leadin','1kHz', '10kHz', '100Hz', 'freqsweep', 'quiet', '3150Hz', '1kHzL', 'swepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV', 'extra_signal','leadout'}; 
signal_names = {'needledrop','leadin','1kHz', '10kHz', '100Hz', 'freqsweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition','leadin2','1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2', '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout2'};
ref_timestamps = [0, 2, 62, 92, 124, 160, 182, 248, 268, 306, 326, 364, 384, 419.5];% this is how many seconds each signal is according to Chris Muth's track listing
lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz

ref_offset = 4.25; % as measured on /020818_A0000B0000/02072019_A0000B000r25-A.wav
ref_timestamps = ref_timestamps + ref_offset;
time = (0:length(reference.dataL)-1)/reference.fs;
ref_transition = 518.25; % as measured on /020818_A0000B0000/02072019_A0000B000r25-A.wav
ref_timestamps2 = ref_timestamps + ref_transition;

%clf(figure(1));
%figure(1);hold on; grid on;
%plot(time, reference.dataL);
%title('Amplitude vs. Time'); 
%xlabel('Time [s]');
%ylabel('Amplitude');


%signal_names = {signal_names, 'transition', signal_names}
ref_signalsL = [];
ref_signalsR = [];
ref_signals = {};

ref_signalsL = reference.dataL(1:floor((ref_timestamps(1))*reference.fs));
ref_signalsR = reference.dataR(1:floor((ref_timestamps(1))*reference.fs));
ref_signals{end + 1} = [ref_signalsL, ref_signalsR];

% first set of signals on the disk
for i = (1:length(ref_timestamps)-1);
    ref_signalsL = reference.dataL(floor(ref_timestamps(i)*reference.fs):floor(ref_timestamps(i+1)*reference.fs));
    ref_signalsR = reference.dataR(floor(ref_timestamps(i)*reference.fs):floor(ref_timestamps(i+1)*reference.fs)); 
    ref_signals{end + 1} = [ref_signalsL, ref_signalsR];
end

% the extended silence section 'transition' between the two sets of signals 
ref_signalsL = reference.dataL(floor(ref_timestamps(end)*reference.fs): ...
                                floor((ref_timestamps(end)+ref_transition)*reference.fs));

ref_signalsR = reference.dataR(floor(ref_timestamps(end)*reference.fs): ...
                                floor((ref_timestamps(end)+ref_transition)*reference.fs));
ref_signals{end + 1} = [ref_signalsL, ref_signalsR];

% second set of signals on the disk
for i = (1:length(ref_timestamps2)-1);
    ref_signalsL = reference.dataL(floor(ref_timestamps2(i)*reference.fs):floor(ref_timestamps2(i+1)*reference.fs));
    ref_signalsR = reference.dataR(floor(ref_timestamps2(i)*reference.fs):floor(ref_timestamps2(i+1)*reference.fs));
    size(ref_signalsL) 
    size(ref_signalsR) 
    ref_signals{end + 1} = [ref_signalsL, ref_signalsR];
end

% get the rest of the file
ref_signalsL = reference.dataL(floor(ref_timestamps2(end)*reference.fs):end);
ref_signalsR = reference.dataR(floor(ref_timestamps2(end)*reference.fs):end);
ref_signals{end + 1} = [ref_signalsL, ref_signalsR];

size(ref_signals)
size(signal_names)

% The loop below is for info only 
for i = (1:length(ref_signals)); 
    disp('next signal');
    size(ref_signals{i});
    length(ref_signals{i})/reference.fs
    %figure(i);
    %plot(ref_signals{i}(:,1));
%    length(ref_signals(i,:,2))
%    length(ref_signals(i,:,1))/reference.fs
%    length(ref_signals(i,:,2))/reference.fs
end


%   newMap = containers.Map(keys, values);
ref_tracks = containers.Map(signal_names, ref_signals);


names = keys(ref_tracks);
datas = values(ref_tracks);

for i = (1:length(ref_tracks));
    %ref_tracks(i)
    disp('next track')
    name = names{i}
    figure(i);
    plot(datas{i});
    title(name);
    length(datas{i})/reference.fs
end

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


