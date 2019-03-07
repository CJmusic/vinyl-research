%{ SCRIPT START BELOW %}
%
path = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.2.wav';
dir_files = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/';
name_files = {'5.1.wav'};  % {['5.2.wav', '5.3.wav', '5.4.wav', '5.5.wav']}; 


[data, time, fs] = audio_load('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.1.wav');

rms(data) %the default matlab rms function works fine

%[data_fft, freq] = audio_spectrum(data, fs, 2^16, 2^16);
%figure(1);
%audio_plotspectrum(freq,data_fft);

figure(2)
%[data_psd, freq_psd] = audio_psd(data, 2^16, fs); 
%audio_plotspectrum(freq_psd, data_psd)
nfft = 2^16
data_L = data(:,1);
%pwelch(data_L, 2^16, fs)
[Pxx,f]=pwelch(data,hanning(nfft,'periodic'),nfft/2,nfft,fs,'');
audio_plotspectrum(f, Pxx, 'PSD')


%    SCRIPT END 

%~~~~~~~~~~~~~~~~~~~~~READ ME~~~~~~~~~~~~~~~~~~~~~~~~~%

%{ This script contains all the functions needed for anaylsis 
%
%
%%} 

%%% HOW TO HANDLE AUDIO ARRAYS in the form [data, time, fs] 


%data = audio_load[1];
%time = audio_load[2];
%fs = audio_load[3];


%data[[LEFT,RIGHT]
%         [LEFT,RIGHT]
%            ...    
%         [LEFT,RIGHT]] 
%left_channel = data[:,1];
%right_channel = data[:,2];

%time = [0.0 sec , 1/fs sec, 2/fs sec, ... , N];


%~~~~~~~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~%

%{ LOAD FILES %}
%This function takes a path and returns the data and time arrays, along with the sample rate
function [data, time, fs] = audio_load(path)
    [data, fs] = audioread(path);
    time = (0:length(data)-1)/fs;
    %data = data.';
end %audio_load

function [data_fft, freq_fft] = audio_spectrum(data, fs, start_sam, n_sam)
    disp('inside audio_spectrum')
    freq_fft = fs*(0:(n_sam/2))/n_sam;
    data_fft = fft(data(start_sam:start_sam+n_sam, :))/n_sam;
    data_fft = data_fft(1:size(data_fft)/2+1);
    disp('finished audio_spectrum')
end %audio_spectrum

function audio_plotspectrum(freq, data_fft, title_string) 
    plot(freq, 20.0*log10(data_fft))  
    set(gca, 'XScale', 'log');
    title(title_string)
    xlabel('Frequency (Hz)')
    ylabel('Level (dB)')  

end

function [data_psd, freq_psd] = audio_psd(data, nfft, fs)
    [data_psd, freq_psd] = pwelch(data , nfft, fs, hanning(nfft), nfft/2);
end


%~ Loop through all the files and line them up with the reference

%{ LINEUP %}
% This function needs to take two arrays, and return two arrays that are lined up with one another. 
% It should: 
%       -  
%
function [cdata, ctime] = lineup(data)
for i = (1:length(name_files));
 
    strcat(dir_files,name_files{i})
    [data_file, fs_file] = audioread(strcat(dir_files,name_files{i}));
    
    data_file = data_file(:,1);
     
    time_file = (0:length(data_file)-1)/fs_file; %Not sure if important, but here I'm using the fs
   
    [acor,lag] = xcorr(data_file,data);
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

function groove_plot(data)
    clf(figure(1))
    clf(figure(2))
    AUDIO_FILES = {};
    
    %~~~Segmenting Parameters~~~%
    time = time;
    rotation_speed = 33.33333;%45;
    T = 60/rotation_speed; %this is the length of one groove segment
    n_sam = round(T*fs)
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
        num_segs = (floor(length(data)/fs/T))
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
end %function plotting
%~~~~~~~~~~~~~~~~~~FUNCTIONS END~~~~~~~~~~~~~~~~~~~~~~%


