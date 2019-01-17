% Record Surface Noise analysis
% using xcorrelation between groove signals
% John Vanderkooy
% Nov. 2018
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
%filename='5.1gaineddeclick.wav';
filename='5.2gaineddeclick.wav';
%filename='mono_noise_20s.wav';
%filename='5wav33type4leadindeclicked.wav';
%filename='5wav_4rev.wav';
%[rev4,fs]=audioread(filename); 
%filename='5wav_4rev_declick.wav';
%filename='5wav_4rev_declick2.wav';
[rev4_declick,fs]=audioread(filename); 
Nt=length(rev4_declick)
Nt_declick=length(rev4_declick) % should be the same...
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%------------------------------------------------
% figure(10)
% plot(t,rev4(:,1),'b') 
% grid on;hold on;
% plot(t,rev4(:,2),'r')
% axis([0 Nt/fs -1 1])
% xlabel('time[s]')
% legend('left','right')
% title('Original')
%-----------------------declicked version------------
figure(12)
plot(t,rev4_declick(:,1),'b') 
grid on;hold on;
plot(t,rev4_declick(:,2),'r')
axis([0 Nt/fs -1 1])
xlabel('time[s]')
legend('left','right')
title('Declicked')
%------------------generate separate groove results-----------
% we'll make n_rev files
n_per_rev=round(fs*60/33.33333)
% n_per_rev=130000
n_rev=floor(Nt/n_per_rev) % number of grooves

% choose the data to plot
for ng=1:n_rev % groove number
    rev4ng(:,:,ng)=rev4_declick(1+(ng-1)*n_per_rev:ng*n_per_rev,:);
end
%---------time variable in plots------------
tg=t(1:n_per_rev);
minplot=-1.0;
maxplot=+1.0;
%-----------------plot signal for each groove-----------------------------
figure(14)
plot(tg,rev4ng(:,1,1),'b') 
grid on;hold on;
plot(tg,rev4ng(:,1,2),'r') 
plot(tg,rev4ng(:,1,3),'g') 
plot(tg,rev4ng(:,1,4),'k') 
axis([0 n_per_rev/fs minplot maxplot])
xlabel('time[s]')
legend('groove1','groove2','groove3','groove4')
title('Declicked left grooves')

%-----------------plot signal for each groove-----------------------------
figure(16)
plot(tg,rev4ng(:,2,1),'b') 
grid on;hold on;
plot(tg,rev4ng(:,2,2),'r'); % variables are sample#, left/right, groove#
plot(tg,rev4ng(:,2,3),'g') 
plot(tg,rev4ng(:,2,4),'k') 
axis([0 n_per_rev/fs minplot maxplot])
xlabel('time[s]')
legend('groove1','groove2','groove3','groove4')
title('Declicked right grooves')
%-----------------plot signal for each groove-----------------------------
figure(18)
plot(tg,rev4ng(:,1,1),'b') 
grid on;hold on;
plot(tg,rev4ng(:,2,1),'r') 
axis([0 n_per_rev/fs minplot maxplot])
xlabel('time[s]')
legend('groove1 left','groove1 right')
title('Declicked stereo groove')

%---------------cross correlation properties-----------
ex1=rev4ng(:,1,1); % groove1 left
ex2=rev4ng(:,1,2); % groove2 left
wy1=rev4ng(:,2,1); % groove1 right
wy2=rev4ng(:,2,2); % groove2 right

Cxy=xcorr(ex1,wy1);
Cxx=xcorr(ex1);
Cyy=xcorr(wy1);
RHOxy=Cxy./(sqrt(Cxx(n_per_rev+1)).*Cyy(n_per_rev+1));

tau=t(1:2*n_per_rev-1);
figure(20)
plot(tau-tau(n_per_rev+1),RHOxy,'b') 
grid on
%axis([0 n_per_rev/fs minplot maxplot])
xlabel('time[s]')
legend('correlation coefficient')
title('Declicked stereo groove')
%----------------------coherence-----------------
nfft=8192;
fc=([0:nfft/2])*fs/nfft;

% rev4ng(:,1,1)=rand(n_per_rev,1);%testing coherence
% rev4ng(:,1,4)=rand(n_per_rev,1);

C=mscohere(rev4ng(:,1,1),rev4ng(:,1,4),hanning(nfft),[],nfft,fs);
%C=mscohere(rev4ng(:,1,1),rev4ng(:,2,1),[],[],nfft,fs);
Cs=boxsmooth(C,0);
figure(25)
plot(fc,C,'b')
axis([0,fs/2,0,1.1])
grid on;
xlabel('frequency [Hz]')
title('Coherence L-R of groove 1')

%C=mscohere(rev4ng(:,1,1),rev4ng(:,1,2),[],[],nfft,fs);
C=mscohere(rev4ng(:,1,1),rev4ng(:,1,2),hanning(nfft),[],nfft,fs);
Cs=boxsmooth(C,0);
figure(26)
plot(fc,C,'b')
axis([0,fs/2,0,1.1])
grid on;
xlabel('frequency [Hz]')
title('Coherence hann L groove1 & L groove2')

C=mscohere(rev4ng(:,1,1),rev4ng(:,2,1),hanning(nfft),[],nfft,fs);
Cs=boxsmooth(C,0);
figure(27)
plot(fc,C,'b')
axis([0,fs/2,0,1.1])
grid on;
xlabel('frequency [Hz]')
title('Coherence L groove1 & L groove4')

C=mscohere(rev4ng(:,2,1),rev4ng(:,2,2),hanning(nfft),[],nfft,fs);
Cs=boxsmooth(C,0);
figure(28)
plot(fc,C,'b')
axis([0,fs/2,0,1.1])
grid on;
xlabel('frequency [Hz]')
title('Coherence R groove1 & R groove2')

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
