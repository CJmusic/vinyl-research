clear all; clc;close all;
disp('----------------start of program--------------------')
set(0,'DefaultLineLinewidth',1.5)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('Volumes/AUDIOBANK/audio_files/A0000B0000')
addpath('Volumes/AUDIOBANK/audio_files/A0137B0137')


trackname = 'sweep';
main_folder = '/Volumes/AUDIOBANK/audio_files/A0000B0000/';
records = {'031418_A0000B0000r028a1558.066','031418_A0000B0000r029a1554.214'};

ref_filename = strcat('/Volumes/AUDIOBANK/audio_files/A0000B0000/031418_A0000B0000r027a1553.770/',trackname,'.wav')


[ref, fs] = audioread(ref_filename);

ts = 10.0;
tf = 20.0;

length(ref)/fs
ref = ref(ts*fs:tf*fs,:);
time = (0:length(ref)-1)/fs;

figure(2)
subplot(2,1,1)
plot(time, ref(:,1))
grid on; hold on;
subplot(2,1,2)
plot(time, ref(:,2))
grid on; hold on;

for i = 1:length(records);
    disp(strcat('i.....', num2str(i)))
    filename = strcat(main_folder,records{i},'/', trackname,'.wav');

    [data, fs] = audioread(filename);

    [acor_L,lags_L] = xcorr(ref(:,1),data(:,1));

    figure(1);
    plot(lags_L/fs,acor_L);
    hold on; grid on;
    [M_L,I_L] = max(abs(acor_L));
    lagdiff_L = lags_L(I_L);    
    lagdiff = lagdiff_L;
    lagdiff = lagdiff/fs;
    lagdiff = 0;

    data = data(floor(ts*fs) + floor(lagdiff*fs):floor(tf*fs) + floor(lagdiff*fs),:);
    
    time = (0:length(data)-1)/fs;

    figure(2)
    subplot(2,1,1)
    plot(time, data(:,1))
    grid on; hold on;
    title('left channel')
    ylim([-1,1])

    subplot(2,1,2)
    plot(time, data(:,2))
    grid on; hold on;
    title('right channel')
    ylim([-1,1])

    nfft=2^14;
    window=hanning(nfft,'periodic');
    [Cc(:,1),fc]=mscohere(ref(:,1),data(:,1),window,nfft/2,nfft,fs);
    [Cc(:,2),fc]=mscohere(ref(:,2),data(:,2),window,nfft/2,nfft,fs);
    [CcLR,fc]=mscohere(data(:,1),data(:,2),window,nfft/2,nfft,fs);


    figure(3)
    subplot(2,1,1)
    plot(Cc(:,1))
    hold on; 
    subplot(2,1,2)
    plot(Cc(:,2))
    hold on; 
    subplot(2,1,1)
    title('left channel')
    xlabel('frequency [Hz]')
    legend('1','2','3','4','5')
    set(gca,'Xscale','log');
    ylim([0,1])
    grid on;
    subplot(2,1,2)
    title('right channel')
    xlabel('frequency [Hz]')
    legend('1','2','3','4','5')
    set(gca,'Xscale','log');
    ylim([0,1])
    grid on;
    % saveas(figure(1), '/Users/cz/Code/vinyl-research/matlab_code/Correlation/plots/coherence.png');





end


