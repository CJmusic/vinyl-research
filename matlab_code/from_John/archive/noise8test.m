% Noise analysis 8 tests
% John Vanderkooy
% Dec 28. 2018
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
filename='EightNoiseTestsDual&Marantz.wav';
[all_eight,fs]=audioread(filename); 
Nt=length(all_eight);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector*******check this
%---------------------all 8 series tests---------------------------
figure(10)
plot(t,all_eight(:,1),'b') 
grid on;hold on;
plot(t,all_eight(:,2),'r')
axis([0 Nt/fs -1 1])
xlabel('time[s]')
legend('left','right')
title('All 8 serial tests')
%------------------generate separate test results-----------
% we'll make 8 files
n_each=20*fs;% 20 seconds for each file
% choose the data to analyse
for ng=1:8 % test number
    test(:,:,ng)=all_eight(1+(ng-1)*n_each:ng*n_each,:);
end
%-------------------time variable in each test------------
tg=t(1:n_each);
minplot=-1.0;
maxplot=+1.0;
%-----------------plot signal for each test-----------------------------
figure(20)
plot(tg,test(:,1,1),'b');% variables are sample#, left/right, test#
grid on;hold on;
plot(tg,test(:,1,2),'r'); 
plot(tg,test(:,1,3),'g'); 
plot(tg,test(:,1,4),'k'); 
plot(tg,test(:,1,5),'b'); 
plot(tg,test(:,1,6),'r'); 
plot(tg,test(:,1,7),'g'); 
plot(tg,test(:,1,8),'k'); 
axis([0 n_each/fs minplot maxplot])
xlabel('time[s]')
%legend('groove1','groove2','groove3','groove4')
title('tests')
%---------------get active portion of each test----------------
ns=floor(11.6*fs);
nf=floor(18.8*fs);
nfft=16384;
for k=1:8;
    d(:,k)=test(ns:nf,1,k);
    [Pxx(:,k) f]=psd(d(:,k),nfft,fs,hanning(nfft),nfft/2);
end
%-----------------plot psd for each test-----------------------------
figure(30)
semilogx(f,10*log10(Pxx(:,1)),'b');
grid on;hold on;
semilogx(f,10*log10(Pxx(:,2)),'r');
semilogx(f,10*log10(Pxx(:,3)),'g');
semilogx(f,10*log10(Pxx(:,4)),'k');
semilogx(f,10*log10(Pxx(:,5)),'b');
semilogx(f,10*log10(Pxx(:,6)),'r');
semilogx(f,10*log10(Pxx(:,7)),'g');
semilogx(f,10*log10(Pxx(:,8)),'k');
%axis([0 n_per_rev/fs minplot maxplot])
xlabel('freq[Hz]');
ylabel('PSD');
legend('shortTS','openTS','phono short','phono open','no table','turntable','leadout','music');
title('PSD of 8 signals')



% %-------------------spectral plots-----------------------
% f=([1:N/2+1]-1)*fs/N;
% f(1)=1e-10;% prevents log plot warning in Octave
% SIL=fft(sil(:,2));
% SILNR=fft(silnr(:,2));
% 
% figure(100);
% semilogx(f,20*log10(abs(SIL(1:N/2+1))),'b');
% grid on;hold on;
% semilogx(f,20*log10(abs(SILNR(1:N/2+1))),'r');
% axis([1 fs/2 -60 40])
% legend('silent track','system noise')
% xlabel('frequency [Hz]');ylabel('dB')
% title(['Spectra']);
% %---------------------------------------
% SILm=pwroctsmooth(SIL,0.1);
% SILNRm=pwroctsmooth(SILNR,0.1);
% 
% figure(200);
% semilogx(f,20*log10(abs(SILm(1:N/2+1))),'b');
% grid on;hold on
% semilogx(f,20*log10(abs(SILNRm(1:N/2+1))),'r');
% axis([1 fs/2 -60 40])
% legend('silent track','system noise')
% xlabel('frequency [Hz]');ylabel('dB')
% title(['Smoothed Spectra']);
% 
% disp('-------------------finished--------------------') 
