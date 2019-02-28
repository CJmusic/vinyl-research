%{ SCRIPT START BELOW %}
%
path_ref = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.2.wav';
dir_files = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/';
name_files = {'5.1.wav'};  % {['5.2.wav', '5.3.wav', '5.4.wav', '5.5.wav']}; 


load_audio('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.1.wav');


%{FUNCTIONS ARE ALL DEFINED BELOW}%


%{ This script contains all the functions needed for anaylsis 
%
%
%%} 

%%% HOW TO HANDLE AUDIO ARRAYS in the form [data_ref, time_ref, fs_ref] 


%data_ref = load_audio[1];
%time_ref = load_audio[2];
%fs_ref = load_audio[3];


%data_ref[[LEFT,RIGHT]
%         [LEFT,RIGHT]
%            ...    
%         [LEFT,RIGHT]] 
%left_channel = data_ref[:,1];
%right_channel = data_ref[:,2];

%time_ref = [0.0 sec , 1/fs sec, 2/fs sec, ... , N];


%~~~~~~~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~%

%{ LOAD FILES %}
%This function takes a path and returns the data and time arrays, along with the sample rate
function [data_ref, time_ref, fs_ref] = load_audio(path_ref)
    [data_ref, fs_ref] = audioread(path_ref);
    data_ref = data_ref(:,1);
    %data_ref = data_ref(1:5*fs_ref);
    size(data_ref)
    time_ref = (0:length(data_ref)-1)/fs_ref;
    AUDIO_FILES{1} = data_ref;
    size(AUDIO_FILES)
    size(name_files)
    name_files
    figure(1); hold on;
    plot(time_ref,data_ref,'g')
end
%~ Loop through all the files and line them up with the reference


%{ LINEUP %}
function [data] = lineup(data)
for i = (1:length(name_files));
 
    strcat(dir_files,name_files{i})
    [data_file, fs_file] = audioread(strcat(dir_files,name_files{i}));
    
    data_file = data_file(:,1);
     
    time_file = (0:length(data_file)-1)/fs_file; %Not sure if important, but here I'm using the fs
    
    [acor,lag] = xcorr(data_file,data_ref);
    [~,I] = max(abs(acor));
    lagDiff = lag(I)
    timeDiff = lagDiff/fs_file
    cdata_file = data_file(lagDiff+1:end);
    ctime_file = (0:length(cdata_file)-1)/fs_file;

    % plot(ctime_file, cdata_file)
    % title('Original Audio Lined up')
    % xlabel('Time (s)')
    % AUDIO_FILES{i+1} = cdata_file; 
end
end

%{ SPECTRUM PLOT %}


%{ GROOVE PLOT %}


clf(figure(1))
clf(figure(2))
AUDIO_FILES = {};

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
name_files= [ {'5.2.wav'},name_files];

figure(2);hold on; legend;
for i=(1:length(AUDIO_FILES))
    data = AUDIO_FILES{i};
    num_segs = (floor(length(data)/fs_ref/T))
    for ng = 1:num_segs
        %seg_array(:,:,ng) = data(1+(ng-1)*n_sam:ng*n_sam,:);
        data_seg = data(1+(ng-1)*n_sam:ng*n_sam,:);
        [s,f,t] = spectrogram(data_seg,winSize,overlap,fftsize,fs_file,'yaxis');
        lo_freq = sum(s(1:5));
        hi_freq = sum(s(6:end));
        %if abs(hi_freq/lo_freq) < 2.0; %%This is where the DETECT SIGNAL algorithm will go once its active
        if ng > 1 &  ng < 4; %for now just plot grooves 2 & 3  
            plot(time_seg,data_seg,  'DisplayName', [ name_files{i},'groove', num2str(ng)]);
        end
    end    
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%~~~~~~~~~~~~~~~~~~FUNCTIONS END~~~~~~~~~~~~~~~~~~~~~~%


