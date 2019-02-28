

path_digital_files = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/00_digital_files/'; 

sweep_file = 'sweep16kHz.wav';
tone100_file = 'tone100.wav';
tone1000_file = 'tone1000.wav'; 
tone3150_file = 'tone3150.wav';
tone10000_file = 'tone10000.wav';

%%{
%
%TRACK LISTINGS ON TEST RECORDS 
%00:02 : 1KHz@7cm-s lateral
%01:02 : 10kHz @ -20dB
%01:32 : 100Hz
%02:04 : frequency sweep 20Hz-16kHz 
%02:40 : quiet groove
%03:02 : 3150 wow & flutter
%04:08 : 1kHz left
%04:28 : frequency sweep left
%05:06 : 1kHz right
%05:26 : frequency sweep right
%06:04 : 1kHz vertical
%06:24 : frequency sweep vertical
%
%%}

%name_files = {[tone1000_file, tone10000_file, tone100_file, sweep_file %{, quiet groove %}, tone3150_file, tone1000_file %{left%}, sweep_file %{left%}, tone1000_file %{right%}, sweep_file %{right%}, tone1000_file %{vertical%}, sweep_file %{vertical%}]};

name_files = {tone1000_file, tone10000_file, tone100_file, sweep_file  tone3150_file, tone1000_file  sweep_file, tone1000_file, sweep_file, tone1000_file, sweep_file};

lagDiff_lists = [0];




%now reference is recording

recorded_file = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/A0000B0000/02072019_A0000B000r25-A.wav'; %path to the recording of the record

[data_rec, fs_rec] = audioread(recorded_file); 

%data_rec = data_rec(:,1); %take only the first channel
%clf(figure(1))
%figure(1) 
%hold on; grid on; 
%plot(time_rec,data_rec, 'g')

[kHztone, fs_file] = audioread(strcat(path_digital_files,name_files{1}));


kHztone = 0.18*kHztone(1:5*fs_file,1);
data_rec = data_rec(6*fs_rec:10*fs_rec,1); %take only the first channel   

time_rec = (0:length(data_rec)-1)/fs_rec;

[acor,lag] = xcorr(kHztone,data_rec);
[~,I] = max(abs(acor));
lagDiff = lag(I)


cdata_file = kHztone(lagDiff+1:end);
ctime_file = (lagDiff:lagDiff + length(cdata_file)-1)/fs_file;


clf(figure(1))
figure(1); hold on; 
plot(time_rec, data_rec) 
plot(ctime_file, cdata_file, 'g')
title('Original Audio Lined up')
xlabel('Time (s)') 




%time_file = (0:length(data_file)-1)/fs_file; %Not sure if important, but here I'm using the fs    
%size(data_rec) 
%size(kHztone)
%[acor,lag] = xcorr(kHztone,data_rec);
%[~,I] = max(abs(acor));
%lagDiff = lag(I)
%timeDiff = lagDiff/fs_file
%%cdata_file = data_file(-lagDiff+1:end);
%%ctime_file = (lagDiff:lagdiff + length(cdata_file)-1)/fs_file;
%lagDiff_lists = [lagDiff_lists, lagDiff]; %append to the list of all lag diffs to keep track of when each signal starts
%plot(ctime_file, cdata_file)
%title('Original Audio Lined up')
%xlabel('Time (s)') 
%AUDIO_FILES{i+1} = cdata_file; 



%
%for i = (1:length(name_files));
%
%    start_offset = lagDiff_lists(end);%[length(lagDiff_lists)];
%    % look at the last entry in the lagDiff, maybe replace this with the ending sample
%     
%    strcat(path_digital_files,name_files{i})
%    [data_file, fs_file] = audioread(strcat(path_digital_files,name_files{i}));
%    
%    data_file = data_file(:,1);
%     
%    time_file = (0:length(data_file)-1)/fs_file; %Not sure if important, but here I'm using the fs
%    
%    [acor,lag] = xcorr(data_file,data_rec);
%    [~,i] = max(abs(acor));
%    lagdiff = lag(i)
%    timediff = lagdiff/fs_file
%    cdata_file = data_file(-lagDiff+1:end);
%    ctime_file = (lagDiff:lagdiff + length(cdata_file)-1)/fs_file;
%    lagDiff_lists = [lagDiff_lists, lagDiff]; %append to the list of all lag diffs to keep track of when each signal starts
%    plot(ctime_file, cdata_file)
%    title('Original Audio Lined up')
%    xlabel('Time (s)')
%    AUDIO_FILES{i+1} = cdata_file; 
%end
%
%
%
