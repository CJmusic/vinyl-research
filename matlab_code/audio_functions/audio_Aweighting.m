% %% A-weighting digital filter
% clear all;close all;clc
% set(0,'DefaultLineLineWidth',1.0);
% set(0,'DefaultAxesFontSize',12);
% set(0,'DefaultAxesFontWeight','bold')
% set(0,'DefaultAxesLineWidth',1.5)

% addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')

% % dBfs = 0;%-34.0%-6.0

% %~~~~~ GENERATE LOG SWEEP ~~~~~~~~%

% % % sig_frac=0.5;% fraction of full scale
% % % sig_frac = 10^(dBfs/20.0)
% % sig_frac = 1.0;
% % fs=96000;% must ensure that Windows settings are the same!
% % N=2^21;%19;% make this larger if there is insufficient time clearance
% % % S_dac=-1.59;% Focusrite 2i2 Yes, it inverts its monitor output!
% % % S_adc=+1.16;% Focusrite 2i2 line input gain @ 12:00 o'clock
% % Npad=0;%N/4;% this is the total zeropad, added to end of play file
% % Ns=N-Npad;% most of array is used for sweep

% % t=linspace(0,(N-1)/fs,N)';% column vector
% % ts=linspace(0,(Ns-1)/fs,Ns)';% to calculate sweep

% % %-----------------------calculate logsweep-----------------------------%
% % f_start=10;% beginning of turnon half-Hann
% % f1=20;% end of turnon half-hann
% % %f2=0.91*fs/2;% beginning of turnoff half-Hann
% % f2=16000%20000;
% % f_stop=fs/2;% end of turnoff half-Hann
% % Ts=Ns/fs;% sweep duration.  This is N-Npad samples in length.
% % Ls=Ts/log(f_stop/f_start);% time for frequency to increase by factor e
% % sweep=zeros(N,1);% initialize
% % sweep(1:Ns)=sin(2*pi*f_start*Ls*(exp(ts/Ls)-1));% logsweep
% % %-----------------------------------------------------------------------%

% % %------------------tapered sweep window------------------------------
% % indexf1=round(fs*Ls*log(f1/f_start))+1;% end of starting taper
% % indexf2=round(fs*Ls*log(f2/f_start))+1;% beginning of ending taper
% % windo=ones(N,1);
% % windo(1:indexf1)=0.5*(1-cos(pi*(1:indexf1)/indexf1));% pre-taper
% % windo(indexf2+1:Ns)=0.5*(1+cos(pi*(1:Ns-indexf2)/(Ns-indexf2)))';% post-taper
% % windo(Ns+1:N)=0; % zeropad end of sweep
% % windosweep=windo.*sweep;% tapered at each end for output to DAC
% % %-----------------------------------------------------------------------%

% % y=sig_frac*[windosweep windosweep];% in phase for vinyl pressing
% %~~~~~ GENERATE LOG SWEEP ENDS ~~~~~~~~%
% fs = 96000;
% time = (1:1/fs:10);
% size(time)
% data(:,1) = sin(2*pi*1000*time);
% data(:,2) = sin(2*pi*1000*time);

% fs = 96000;
% figure(1);
% plot(data)
% grid on;

% [data_fft, freq_fft] = audio_spectrum(data, fs, 1, length(data));
% figure(2)
% audio_plotspectrum(freq_fft, data_fft, 'pre filter')


% % data = data(1:10*fs);
% data_A = audio_AweightingTest(data);

% [dataA_fft, freqA_fft] = audio_spectrum(data_A, fs, 1, length(data_A));
% figure(3)
% audio_plotspectrum(freqA_fft, dataA_fft, 'post filter')

% tracks = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r029a.wav');
% tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/-01a.wav');
% data = tracks('1kHz');

% fs = 96000;
% figure(4);
% plot(data)
% grid on;

% [data_fft, freq_fft] = audio_spectrum(data, fs, 1, length(data));
% figure(5)
% audio_plotspectrum(freq_fft, data_fft, 'pre filter')


% % data = data(1:10*fs);
% data_A = audio_AweightingTest(data);

% [dataA_fft, freqA_fft] = audio_spectrum(data_A, fs, 1, length(data_A));
% figure(6)
% audio_plotspectrum(freqA_fft, dataA_fft, 'post filter')

function data_A = audio_Aweighting(data)
% function data_A = audio_AweightingTest(data)
    fs = 96000;
    f1 = 20.598997; 
    f2 = 107.65265;
    f3 = 737.86223;
    f4 = 12194.217;%for infinite fs
    f4 = 14100;%for finite fs, makes magnitude -9.5dB at 20kHz
    A1000 = 1.9997;
    % A1000 = 0.0001;
    NUM = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
    DEN = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
    DEN = conv(conv(DEN,[1 2*pi*f3]),[1 2*pi*f2]);
    % % Bilinear transformation of analog design to get the digital filter. 
    [b,a] = bilinear(NUM,DEN,fs);
    % time = (0:length(data)-1)/fs;
    data_A = filter(b, a, data);
    % data_A = filtfilt(b, a, data);
end