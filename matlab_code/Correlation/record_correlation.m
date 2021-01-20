% Record Surface Noise analysis
% using correlation & coherence between different recordings same groove
% John Vanderkooy
% Feb. 2019
%  
clear all; clc;close all;
disp('----------------start of program--------------------')
set(0,'DefaultLineLinewidth',1.5)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
%
try
    pkg load signal %for Octave
catch
end

addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')


trackname = 'transition'

reftracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/SameRecordTest2/r300a1550.844.wav');
% folder = '/Volumes/AUDIOBANK/audio_files/SameRecordTest2/';
fs = 96000; 

[b,a]=butter(2,2*200/fs,'high');% remove LF arm resonance


reftracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/coherencetest/A0000B0000/031418_A0000B0000r027a1553.770.wav')
folder = '/Volumes/AUDIOBANK/audio_files/coherencetest/A0000B0000';

files = dir(fullfile(folder,'*.wav'))

for i = (1:length(files)) %%loop through records
    filename = files(i).name;
    tracks = SeperateTracks(strcat(files(i).folder,'/',filename));
    data = tracks(trackname);
    ref = reftracks(trackname);

    if length(data)<length(ref)
        ref = ref(1:length(data),:);
    end

    if length(ref)<length(data)
        data = data(1:length(ref),:);
    end

    [acor_L,lags_L2] = xcorr(ref(:,1),data(:,1));
    [M_L,I_L] = max(abs(acor_L));
    lagdiff_L2 = lags_L2(I_L);
    lagdiff2 = lagdiff_L2;
    data = circshift(data,lagdiff2);

    % maxlag=2^16;
    % dataf=filter(b,a,data);
    % XCORR(:,1) = xcorr(reff(:,1),dataf(:,1),maxlag);
    % XCORR(:,2) = xcorr(reff(:,2),dataf(:,2),maxlag);
    % time = (0:(length(data)-1)/fs);
    % time = t;
    nfft=2^14;
    window=hanning(nfft,'periodic');
    [PSD(:,1),f]=pwelch(data(:,1),window,nfft/2,nfft,fs,'psd');
    [PSD(:,2),f]=pwelch(data(:,2),window,nfft/2,nfft,fs,'psd');
    
    [Cc(:,1),fc]=mscohere(ref(:,1),data(:,1),window,nfft/2,nfft,fs);
    [Cc(:,2),fc]=mscohere(ref(:,2),data(:,2),window,nfft/2,nfft,fs);

    [data_fft_L, freq] = audio_spectrum(data(:,1), fs, 1, length(data));
    [data_fft_R, freq] = audio_spectrum(data(:,2), fs, 1, length(data));
    
    figure(1)
    subplot(2,1,1)
    plot(Cc(:,1))
    hold on; 
    subplot(2,1,2)
    plot(Cc(:,2))
    hold on; 

    figure(2)
    subplot(2,1,1)
    plot(PSD(:,1))
    hold on; 
    subplot(2,1,2)
    plot(PSD(:,2))
    hold on; 

    time=(0:length(data(:,1))-1)/fs;%column vector
    figure(3)
    subplot(2,1,1)
    plot(time,data(:,1))
    hold on;
    subplot(2,1,2)
    time=(0:length(data(:,2))-1)/fs;%column vector
    plot(time,data(:,2))
    hold on;

    figure(4)
    subplot(2,1,1)
    plot(freq, 20.0*log10(abs(data_fft_L))) 
    grid on; hold on;
    set(gca, 'XScale', 'log');
    xlabel('Frequency (Hz)')
    ylabel('Level (dB)')  
    subplot(2,1,2)
    plot(freq, 20.0*log10(abs(data_fft_R))) 
    grid on ; hold on;
    set(gca, 'XScale', 'log');
    xlabel('Frequency (Hz)')
    ylabel('Level (dB)')  
    % size(time)
    % size(data)

    % figure(2)
    % plot(time, data)
    % hold on;



end


figure(1)
grid on;
% set(gca,'Xscale','log');
% title('Coherence of multiple recordings of the same record')
subplot(2,1,1)
title('left channel')
xlabel('frequency [Hz]')
legend('1','2','3','4','5')
set(gca,'Xscale','log');
grid on;
subplot(2,1,2)
title('right channel')
xlabel('frequency [Hz]')
legend('1','2','3','4','5')
set(gca,'Xscale','log');
grid on;
saveas(figure(1), '/Users/cz/Code/vinyl-research/matlab_code/Correlation/plots/coherence.png');


figure(2)
grid on;
% set(gca,'Xscale','log');
% title('Coherence of multiple recordings of the same record')
subplot(2,1,1)
title('left channel')
xlabel('frequency [Hz]')
legend('1','2','3','4','5')
set(gca,'Xscale','log');
grid on;
subplot(2,1,2)
title('right channel')
xlabel('frequency [Hz]')
legend('1','2','3','4','5')
set(gca,'Xscale','log');
grid on;
saveas(figure(2), '/Users/cz/Code/vinyl-research/matlab_code/Correlation/plots/PSD.png');

 

figure(3)
subplot(2,1,1)
title('left channel')
legend('1','2','3','4','5')
ylim([-1,1])
grid on;
subplot(2,1,2)
title('right channel')
legend('1','2','3','4','5')
ylim([-1,1])
grid on;
saveas(figure(3), '/Users/cz/Code/vinyl-research/matlab_code/Correlation/plots/signal.png');

figure(4)
subplot(2,1,1)
title('left channel')
legend('1','2','3','4','5')
ylim([-196,0])
grid on;
subplot(2,1,2)
title('right channel')
legend('1','2','3','4','5')
ylim([-196,0])
grid on;
saveas(figure(4), '/Users/cz/Code/vinyl-research/matlab_code/Correlation/plots/spectrum.png');

