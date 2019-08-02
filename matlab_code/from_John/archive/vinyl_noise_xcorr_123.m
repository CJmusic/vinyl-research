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
filename='Bcorrelation_test_1.wav';ts=7.2028;tf=26.5;
[rev4_1,fs]=audioread(filename);
filename='Bcorrelation_test_2.wav';ts2=11.50612;
[rev4_2,fs]=audioread(filename);
filename='Bcorrelation_test_3.wav';ts3=16.33612;
[rev4_3,fs]=audioread(filename);
lr=1;%1=left, 2=right
Nt=length(rev4_1);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration_1: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-------------------plot all data------------
figure(20)
plot(t,rev4_1(:,1),'b') 
grid on;hold on;
% plot(t,rev4_2(:,1),'r')
% plot(t,rev4_3(:,1),'g')
% axis([0 Nt/fs ylim])
% xlabel('time[s]')
% legend('left','right')
% title('all recorded data')
%----------------choose equal length files----------------
% ts and tf should be determined from figure 20
%tf=8.0;% use this to shorten all the runs
ns=round(ts*fs);nf=round(tf*fs);
rev4_1=rev4_1(ns:nf,:);
Nt=length(rev4_1);
ns2=round(ts2*fs);
rev4_2=rev4_2(ns2:ns2+Nt-1,:);
ns3=round(ts3*fs);
rev4_3=rev4_3(ns3:ns3+Nt-1,:);
t=linspace(0,(Nt-1)/fs,Nt)';%column vector

figure(30)
plot(t,rev4_1(:,1),'b') 
grid on;hold on;
plot(t,rev4_2(:,1),'r')
plot(t,rev4_3(:,1),'g')
axis([0 Nt/fs ylim])
xlabel('time[s]')
legend('1','2','3')
title('all recorded data')
%---------------histogram----------------
figure(35)
hist(rev4_1,2048)
grid on;
title('histogram')
%--------------psd of whole segment----------------
nfft=2^13;
[Pxx,f]=pwelch(rev4_1(:,lr),hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
[Pxx2,f]=pwelch(rev4_2(:,lr),hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
[Pxx3,f]=pwelch(rev4_3(:,lr),hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
psdrms=sqrt(sum(Pxx)*fs/nfft)
%---------------smoothing if necessary----------------
% Pxxc=pwroctsmooth(Pxx,0.1);
%-----------------plot psd for each test-----------------------------
figure(40)
semilogx(f,10*log10(Pxx),'b');
grid on;hold on;
semilogx(f,10*log10(Pxx2),'r');
semilogx(f,10*log10(Pxx3),'g');
axis([xlim -140 -40])
xlabel('freq[Hz]')
ylabel('Power Spectrum')
legend('1','2','3','Location','Best');
title('PSD')
%----------------------coherence-----------------
[Cc1(:,lr),fc]=mscohere(rev4_1(:,lr),rev4_2(:,lr),hanning(nfft),nfft/2,nfft,fs);
[Cc2(:,lr),fc]=mscohere(rev4_1(:,lr),rev4_3(:,lr),hanning(nfft),nfft/2,nfft,fs);
%[Cc3(:,lr),fc]=mscohere(rev4_2(:,lr),rev4_3(:,lr),hanning(nfft),nfft/2,nfft,fs);
%---------------plot coherences---------
figure(50)
semilogx(fc,Cc1(:,lr),'r')
grid on;hold on;
semilogx(fc,Cc2(:,lr),'g')
%semilogx(fc,Cc3(:,lr),'b')
axis([0,fs/2,0,1])
legend('blu-red','blu-grn','red-grn')
xlabel('frequency [Hz]')
title('Run Coherences')

figure(60)
plot(fc,Cc1(:,lr),'r')
grid on;hold on;
plot(fc,Cc2(:,lr),'g')
%plot(fc,Cc3(:,lr),'b')
axis([0,0.05*fs/2,0,1])
legend('blu-red','blu-grn','red-grn','location','best')
xlabel('frequency [Hz]')
title('Run Coherences')
disp('-------------------finished--------------------') 
