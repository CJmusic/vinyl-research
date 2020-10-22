% %% CCIR-weighting digital filter
% clear all;close all;clc
% set(0,'DefaultLineLineWidth',1.0);
% set(0,'DefaultAxesFontSize',12);
% set(0,'DefaultAxesFontWeight','bold')
% set(0,'DefaultAxesLineWidth',1.5)

% addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')

% % dBfs = 0;%-34.0%-6.0

% %~~~~~~~ GENERATE LOG SWEEP ~~~~~~~~%

% % sig_frac=0.5;% fraction of full scale
% % sig_frac = 10^(dBfs/20.0)
% sig_frac = 1.0;
% fs=96000;% must ensure that Windows settings are the same!
% N=2^21;%19;% make this larger if there is insufficient time clearance
% % S_dac=-1.59;% Focusrite 2i2 Yes, it inverts its monitor output!
% % S_adc=+1.16;% Focusrite 2i2 line input gain @ 12:00 o'clock
% Npad=0;%N/4;% this is the total zeropad, added to end of play file
% Ns=N-Npad;% most of array is used for sweep

% t=linspace(0,(N-1)/fs,N)';% column vector
% ts=linspace(0,(Ns-1)/fs,Ns)';% to calculate sweep

% %-----------------------calculate logsweep-----------------------------%
% f_start=10;% beginning of turnon half-Hann
% f1=20;% end of turnon half-hann
% %f2=0.91*fs/2;% beginning of turnoff half-Hann
% f2=16000%20000;
% f_stop=fs/2;% end of turnoff half-Hann
% Ts=Ns/fs;% sweep duration.  This is N-Npad samples in length.
% Ls=Ts/log(f_stop/f_start);% time for frequency to increase by factor e
% sweep=zeros(N,1);% initialize
% sweep(1:Ns)=sin(2*pi*f_start*Ls*(exp(ts/Ls)-1));% logsweep
% %-----------------------------------------------------------------------%

% %------------------tapered sweep window------------------------------
% indexf1=round(fs*Ls*log(f1/f_start))+1;% end of starting taper
% indexf2=round(fs*Ls*log(f2/f_start))+1;% beginning of ending taper
% windo=ones(N,1);
% windo(1:indexf1)=0.5*(1-cos(pi*(1:indexf1)/indexf1));% pre-taper
% windo(indexf2+1:Ns)=0.5*(1+cos(pi*(1:Ns-indexf2)/(Ns-indexf2)))';% post-taper
% windo(Ns+1:N)=0; % zeropad end of sweep
% windosweep=windo.*sweep;% tapered at each end for output to DAC
% %-----------------------------------------------------------------------%

% y=sig_frac*[windosweep windosweep];% in phase for vinyl pressing
% %~~~~~ GENERATE LOG SWEEP ENDS ~~~~~~~~%
% fs = 96000;
% time = (1:1/fs:10);
% size(time)
% data(:,1) = sin(2*pi*2000*time);
% data(:,2) = sin(2*pi*2000*time);

% fs = 96000;
% figure(1);
% plot(data)
% grid on;

% [data_fft, freq_fft] = audio_spectrum(data, fs, 1, length(data));
% figure(2)
% audio_plotspectrum(freq_fft, data_fft, 'pre filter')


% % data = data(1:10*fs);
% data_A = audio_CCIRweightingTest(data);

% [dataA_fft, freqA_fft] = audio_spectrum(data_A, fs, 1, length(data_A));
% figure(3)
% audio_plotspectrum(freqA_fft, dataA_fft, 'post filter')

% % tracks = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r029a.wav');
% % tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/-01a.wav');
% % data = tracks('1kHz');

% data = y;

% fs = 96000;
% figure(4);
% plot(data)
% grid on;

% [data_fft, freq_fft] = audio_spectrum(data, fs, 1, length(data));
% figure(5)
% audio_plotspectrum(freq_fft, data_fft, 'pre filter')


% % data = data(1:10*fs);
% data_A = audio_CCIRweightingTest(data);

% [dataA_fft, freqA_fft] = audio_spectrum(data_A, fs, 1, length(data_A));
% figure(6)
% audio_plotspectrum(freqA_fft, dataA_fft, 'post filter')

function data_CCIR = audio_CCIRweighting(data)
% function data_CCIR = audio_CCIRweightingTest(data)
   fs=96000;
   frdc=[0 31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000 25000 30000 fs/2];
   CCIR=[-inf -35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8 -35 -50 -inf];
   Wn=2*frdc/fs;
   CCIRmag=10.^(CCIR/20);
   [b,a]=yulewalk(12,Wn,CCIRmag);%%%%%%%%%%%%%%%%%%%%%
   [d,c]=butter(1,2*750/fs,'high');% this corrects DC-LF with highpass
   fb=conv(b,d);ea=conv(a,c);
   data_CCIR=filter(fb,ea,data);
   % data_CCIR=filtfilt(fb,ea,data);
   %----------------------plot yulewalk freq------------
   % f=[0:N/2]*fs/N;
   % OUTPUT=fft(signal);
   %------------------------------------------------------------------
   
end