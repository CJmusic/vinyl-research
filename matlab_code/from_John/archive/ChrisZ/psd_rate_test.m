% PSD sample rate variation test
% John Vanderkooy
% Feb. 2019
% the tricky thing here is to know what we actually want to calculate
% let us calculate the noise as preceived for the same signal length,
% thus this is like continuous noise heard by a listener, ie, the psd
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
%---------------implementation constants-------------------
N=2^20;% let us use 1/4 of this as the original, thus alphamax=4
fs=44100;%original sampling rate
%---------------------------------------------
sig=rand(N,1)-0.5;
%----------------filter data to give psd some character-----------------
[b,a]=butter(6,[0.1 0.4]);
sigf=filter(b,a,sig);
%sigf=sig;%---------!!!!!!!!!!!!!
rms=sqrt(sum(sigf.^2)/N)
%---------------use and plot 1/4 of original data--------
Nt=N/4;
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-----------------plot original sampled signal---------------
figure(10)
plot(t,sigf(1:Nt),'b') 
grid on;
axis([0 t(Nt) ylim])
xlabel('time[s]')
ylabel('amplitude')
title('original data (N/4)')
%-----------------rate shifted data-------------
P=4;Q=2;alpha=P/Q;% resampling freq ratio
sigshf=resample(sigf,P,Q);
fs2=alpha*fs;
Nsh=floor(alpha*Nt);% new samples with same duration as N/4 original
sigshf=alpha*sigf(1:Nsh);% assumes the data results from a differentiating device
%sigshf=sigf(1:Nsh);% assumes the data is normal, power preserved 
tsh=linspace(0,(Nsh-1)/fs2,Nsh)';
rms=sqrt(sum(sigshf.^2)/Nsh)% should be up by alpha ? but psd ?
%--------------------plot speeded up filtered signal------------
figure(20)
plot(tsh,sigshf,'b') 
grid on;
axis([0 tsh(Nsh) ylim])
xlabel('Time[s]')
title('rate speeded data')
%--------------power spectral density plots--------------------
nfft=1024;
%the matlab psd function may not be correct
[Pxx,f]=pwelch(sigf,hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
psdrms=sqrt(sum(Pxx)*fs/nfft)
[Pxx2,fsh]=pwelch(sigshf,hanning(nfft,'periodic'),nfft/2,nfft,fs2,'psd');
psdrms=sqrt(sum(Pxx2)*fs2/nfft)
%---------------smoothing if necessary----------------
% Pxxc=pwroctsmooth(Pxx,0.2);
% Pxx2c=pwroctsmooth(Pxx2,0.2);
%-----------------plot psd for each test-----------------------------
figure(30)
semilogx(f,10*log10(Pxx),'b');
grid on;hold on;
semilogx(fsh,10*log10(Pxx2),'r');
axis([xlim -100 -40])
xlabel('freq[Hz]');
ylabel('Power Spectrum');
legend('normal psd' ,'rate speeded psd','Location','Best');
title('PSD')
disp('-------------------finished--------------------') 
