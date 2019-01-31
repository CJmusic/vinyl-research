%reference = 'correlation1.wav' 
%filename = 'correlation2.wav'
%
%
%filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.1.wav';
%reference = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.2.wav';
%[data_file, fs] = audioread(filename);
%[data_ref, fs] = audioread(reference);
%
%data_ref = data_ref((length(data_ref)/3):length(data_ref));
%data_file = data_file((length(data_file)/3):length(data_file));
%
%t_ref = (0:length(data_ref)-1)/fs;
%t_file = (0:length(data_file)-1)/fs;
%
%figure(2)
%subplot(2,1,1)
%plot(t_ref,data_ref)
%title('s_1')
%
%subplot(2,1,2)
%plot(t_file,data_file)
%title('s_2')
%xlabel('Time (s)')
%
%[acor,lag] = xcorr(data_file,data_ref);
%
%[~,I] = max(abs(acor));
%lagDiff = lag(I)
%
%timeDiff = lagDiff/fs
%
%%figure(3)
%%plot(lag,acor)
%
%cdata_file = data_file(lagDiff+1:end); %might be negative -lagDiff
%ct_file = (0:length(cdata_file)-1)/fs;
%
%figure(1)
%subplot(2,1,1)
%plot(t_ref,data_ref)
%title('s_1, aligned')
%%xlim([0,1])
%subplot(2,1,2)
%
%plot(ct_file,cdata_file)
%title('s_2')
%xlabel('Time (s)')
%xlim([0,1])


clf(figure(1))
clf(figure(2))

%{
This code uses xcorr to lineup two audio files. The code above is the simplest use case, below it's been specialized
to work with multiple recordings of the same audio file and our groove plotting method.
%}

path_ref = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.2.wav';
dir_files = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/';
name_files = {'5.1.wav'};  % {['5.2.wav', '5.3.wav', '5.4.wav', '5.5.wav']}; 




%path_ref = 'correlation1.wav' 
%dir_files = ''
%name_files = {'correlation2.wav'}

AUDIO_FILES = {};

[data_ref, fs_ref] = audioread(path_ref);
data_ref = data_ref(:,1);
%data_ref = data_ref(1:5*fs_ref);
size(data_ref)
time_ref = (0:length(data_ref)-1)/fs;
AUDIO_FILES{1} = data_ref;
size(AUDIO_FILES)
size(name_files)
name_files
%~ Loop through all the files and line them up with the reference
figure(1); hold on;
plot(time_ref,data_ref, 'g')
title('s_1')



%%Gonna slice out a smaller section of audio to line them up

lineup_ref  = data_ref;%(10*fs_ref:15*fs_ref);
size(lineup_ref)
for i = (1:length(name_files));
 
    strcat(dir_files,name_files{i})
    [data_file, fs_file] = audioread(strcat(dir_files,name_files{i}));
    
    data_file = data_file(:,1);
    %data_file = data_file(1:7*fs_ref);
     
    lineup_file = data_file;%(10*fs_file:15*fs_file); 
    size(lineup_ref)
   
    time_file = (0:length(data_file)-1)/fs_file; %Not sure if important, but here I'm using the fs
    
    [acor,lag] = xcorr(lineup_file,lineup_ref);
    [~,I] = max(abs(acor));
    lagDiff = lag(I)
    timeDiff = lagDiff/fs_file

    cdata_file = data_file(lagDiff+1:end);
    ctime_file = (0:length(cdata_file)-1)/fs_file;
   
    size(data_file)
    size(cdata_file)
    size(data_ref)
   
    plot(ctime_file, cdata_file)
    title('s_2')
    xlabel('Time (s)')
    AUDIO_FILES{i+1} = cdata_file; 
end

%~~~Segmenting Parameters~~~%
time = time_ref;
rotation_speed = 33.33333;%45;
T = 60/rotation_speed; %this is the length of one groove segment
n_sam = round(T*fs_ref)
time_seg = time(1:n_sam);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%~~Detect_Signal Parameters~~~%
winSize= 2^8;%round(fs*dur);
overlap=0;%round(winSize/2);
fftsize=2^8;%winSize;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%~~~~~GROOVE PLOTTING~~~~~~%
seg_array = [];

figure(2);hold on;

for i=(1:length(AUDIO_FILES))
    data = AUDIO_FILES{i};
    num_segs = (floor(length(data)/fs_ref/T))
    for ng = 1:num_segs
        %seg_array(:,:,ng) = data(1+(ng-1)*n_sam:ng*n_sam,:);
        %%IF THERE HI FREQ DONT PLOT IT
        data_seg = data(1+(ng-1)*n_sam:ng*n_sam,:);
        [s,f,t] = spectrogram(data_seg,winSize,overlap,fftsize,fs,'yaxis');
        lo_freq = sum(s(1:5));
        hi_freq = sum(s(6:end));
        if abs(hi_freq/lo_freq) < 1.0;
            plot(time_seg,data_seg);
        end
    end    
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~%
