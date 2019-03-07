
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/260219_noisereferenceinst/');


%addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/');
%audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/'
%
%AUDIO_FILES = {'one.wav','two.wav','three.wav', 'four.wav', 'five.wav'};

addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/');
audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/';

AUDIO_FILES = {'Bcorrelation_test_1.wav','Bcorrelation_test_2.wav','Bcorrelation_test_3.wav'};

%AUDIO_FILES = {[fsinst_shorted, fsinst_shortedgained, gain10reference, gain10recordnoise, recordnoise, reference, system_noise, system_gained] };

path_ref = strcat(audio_dir,AUDIO_FILES{1});
path_ref
[data_ref, fs_ref] = audioread(path_ref);
time_ref = (0:length(data_ref)-1)/fs_ref;

ref_coh = data_ref;%(length(data_ref)/4:length(data_ref)/4+2^20,:);

figure(1);
grid on; hold on;
plot(time_ref,data_ref,'g');

for i = (1:length(AUDIO_FILES)); 
    strcat(audio_dir,AUDIO_FILES{i})
    [data, time, fs] = audio_load(strcat(audio_dir,AUDIO_FILES{i}));
    [cdata, ctime] = audio_lineup(data, fs, time, data_ref);

    file_coh = cdata;%(length(cdata)/4:length(cdata)/4+2^16,:);


    [ amp_coh, freq_coh ] = audio_mscohere(ref_coh, file_coh, fs_ref);
%    nfft=2^14;
%    freq_coh=([0:nfft/2])*fs_ref/nfft;
%    amp_coh = mscohere(ref_coh,file_coh,hanning(nfft),[],nfft,fs_ref);
    
    figure(1);
    plot(ctime,cdata)

    figure(2)
    %semilogx(freq_coh,amp_coh(:,1))
    plot(freq_coh,amp_coh(:,1))
    grid on;hold on;
    axis([0,fs_ref/2,0,1])
    xlabel('frequency [Hz]')
    title('Groove Coherences')
    
end 

%%----- COHERENCE -----%
%nfft=2^14;
%freq_coh=([0:nfft/2])*fs_ref/nfft;
%amp_coh = mscohere(data_ref,data_ref,hanning(nfft),[],nfft,fs_ref);
%semilogx(freq_coh,amp_coh(:,1),'b')
%grid on;hold on;
%axis([0,fs_ref/2,0,1])
%xlabel('frequency [Hz]')
%title('Groove Coherences')
%%---------------------%


%~~~ Loop through all the files and line them up with the reference ~~~%
%for i = (1:length(name_files));
% 
%    strcat(dir_files,name_files{i})
%    [data_file, fs_file] = audioread(strcat(dir_files,name_files{i})); 
%    data_file = data_file(:,1); %Take only the first channel (L) 
%    time_file = (0:length(data_file)-1)/fs_file; %Not sure if important, but here I'm using the fs
%    
%    [acor,lag] = xcorr(data_file,data_ref);
%    [~,I] = max(abs(acor));
%    lagDiff = lag(I)
%    timeDiff = lagDiff/fs_file
%    cdata_file = data_file(lagDiff+1:end);
%    ctime_file = (0:length(cdata_file)-1)/fs_file;
%    plot(ctime_file, cdata_file)
%    title('Original Audio Lined up')
%    xlabel('Time (s)')
%    AUDIO_FILES{i+1} = cdata_file; 
%end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%%~~~Segmenting Parameters~~~%
%time = time_ref;
%rotation_speed = 33.33333;%45;
%T = 60/rotation_speed; %this is the length of one groove segment
%n_sam = round(T*fs_ref)
%time_seg = time(1:n_sam);
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%
%%~~~~~SEGMENTING~~~~~%
%%Remap the data file to a single array parsed out via groove number 
%for ng = 1:num_segs
%    seg_array(:,:,ng) = data(1+(ng-1)*n_sam:ng*n_sam,:);
%end
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%
%
%%~~Detect_Signal Parameters~~~%
%winSize= 2^8;%round(fs*dur);
%overlap=0;%round(winSize/2);
%fftsize=2^8;%winSize;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%
%
%%~~~~~GROOVE, COHERANCE and MSCOHERE PLOTTING~~~~~~%
%seg_array = [];
%name_files= [ {'5.2.wav'},name_files];
%
%figure(2);hold on; legend;
%
%data = data_ref ;% AUDIO_FILES{i};
%num_segs = (floor(length(data)/fs_ref/T))
%data_groove = data(1:n_sam,:);
%
%for ng = 1:num_segs-3
%    %seg_array(:,:,ng) = data(1+(ng-1)*n_sam:ng*n_sam,:);
%    data_seg = data(1+(ng-1)*n_sam:ng*n_sam,:);
%
%    corr_coeff = sum(rev4ng(:,1,1).*rev4ng(:,1,ng),1);
%
%    figure(2);hold on; legend;
%    plot(time_seg,data_seg,  'DisplayName', [ name_files{i},'groove', num2str(ng)]);
%    
%    [cxy, f] = mscohere(data(1:n_sam,:),data_seg,fs_ref,[]);
%
%    figure(3); hold on; legend; 
%    plot(f*fs_ref, cxy);
%    set(gca, 'XScale', 'log');
%end   
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%
%%
%%correlation function from John
%for ng=1:n_rev
%    cf(ng)=sum(rev4ng(:,1,1).*rev4ng(:,1,ng),1);
%end
%figure(13)
%plot(cf)
%xlabel('groove number')
%ylabel('correlation coefficient')
%grid on
