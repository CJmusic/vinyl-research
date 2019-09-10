% Record Surface Noise analysis
% using correlation & coherence between groove signals
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
%filename='Bsilent_cartridge_96_no_w_load.wav';ts=0;tf=13;% no load
filename='Bsilent_cartridge_96_no_w_load.wav';ts=14;tf=26.0;% 47K load
%filename='Bnosig_cartridge_96.wav';ts=0.1;tf=20.0;
%filename='Bsilent_cartridge_96.wav';ts=0.1;tf=20.0;
%filename='ChirpSilentLeadin.wav';ts=0;tf=21;%silent closely spaced grooves
%filename='mono_noise_20s.wav';ts=0;tf=19;
[rev4,fs]=audioread(filename);
lr=1;%1=left, 2=right
Nt=length(rev4);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-------------------plot all data------------
figure(20)
plot(t,rev4(:,1),'b') 
grid on;hold on;
plot(t,rev4(:,2),'r')
axis([0 Nt/fs ylim])
xlabel('time[s]')
legend('left','right')
title('all recorded data')
%----------------choose data to analyze----------------
% ns and nt should be defined with filename above
ns=round(ts*fs)+1;nf=round(tf*fs);
if tf>Nt/fs
    tf=Nt/fs;ts=0;
end
rev4=rev4(ns:nf,:);
Nt=length(rev4)
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%---------------histogram----------------
figure(25)
hist(rev4,2048)
grid on;
title('histogram')
%--------------psd of whole segment----------------
nfft=2^12;
%the matlab psd function may not be correct. pwelch seems OK
[Pxx,f]=pwelch(rev4(:,lr),hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
psdrms=sqrt(sum(Pxx)*fs/nfft)
%---------------smoothing if necessary----------------
% Pxxc=pwroctsmooth(Pxx,0.1);
%-----------------plot psd for each test-----------------------------
figure(26)
semilogx(f,10*log10(Pxx),'b');
grid on;
% hold on;
% semilogx(fsh,10*log10(Pxx2),'r');
axis([xlim -150 -90])
xlabel('freq[Hz]')
ylabel('Power Spectrum')
legend('normal psd')% ,'rate speeded psd','Location','Best');
title('PSD')

%------------------generate separate groove results-----------
n_per_rev=round(fs*60/33.3333333)
n_rev=floor(Nt/n_per_rev) % number of grooves
rev4ng=zeros(n_per_rev,2,n_rev);
for ng=1:n_rev % groove number
    rev4ng(:,:,ng)=rev4(1+(ng-1)*n_per_rev:ng*n_per_rev,:);
end
%---------time variable in plots------------
tg=t(1:n_per_rev);
%-----------------plot signal for each groove-----------------------------
figure(30)
plot(tg,rev4ng(:,1,1),'b') 
grid on;hold on;
for k=2:n_rev
    plot(tg,rev4ng(:,1,k),'r'); % variables are sample#, left/right, groove#
end
axis([0 t(n_per_rev) ylim])
xlabel('time[s]')
title('all left grooves')
%-----------------plot signal for each groove-----------------------------
figure(40)
plot(tg,rev4ng(:,2,1),'b') 
grid on;hold on;
for k=2:n_rev
    plot(tg,rev4ng(:,2,k),'r'); % variables are sample#, left/right, groove#
end
axis([0 n_per_rev/fs ylim])
xlabel('time[s]')
title('all right grooves')
%---------------correlation between grooves--------------
for k=1:n_rev
C(k)=sum(rev4ng(:,lr,1).*rev4ng(:,lr,k));
end
figure(50)
plot(C)
grid on
xlabel('groove #')
axis([1 n_rev 0 C(1)])
title('correlation coefficient')
%-----------------plot signal for each groove-----------------------------
figure(60)
plot(tg,rev4ng(:,1,1),'b') 
grid on;hold on;
plot(tg,rev4ng(:,2,1),'r') 
axis([0 n_per_rev/fs ylim])
xlabel('time[s]')
legend('groove1 left','groove1 right')
title('stereo groove')
%------------------xcorr-----------------
% Xc=xcorr(rev4ng(:,1,1),rev4ng(:,2,1));
% figure(65)
% plot(tg-n_per_rev/fs,Xc)

%---------------cross correlation properties-----------
% ex1=rev4ng(:,1,1); % groove1 left
% ex2=rev4ng(:,1,2); % groove2 left
% wy1=rev4ng(:,2,1); % groove1 right
% wy2=rev4ng(:,2,2); % groove2 right
% 
% Cxy=xcorr(ex1,wy1);
% Cxx=xcorr(ex1);
% Cyy=xcorr(wy1);
% RHOxy=Cxy./(sqrt(Cxx(n_per_rev+1)).*Cyy(n_per_rev+1));
% 
% tau=t(1:2*n_per_rev-1);
% figure(20)
% plot(tau-tau(n_per_rev+1),RHOxy,'b') 
% grid on
% %axis([0 n_per_rev/fs minplot maxplot])
% xlabel('time[s]')
% legend('correlation coefficient')
% title('Declicked stereo groove')
%----------------------coherence-----------------
nfft=2^14;
fc=([0:nfft/2])*fs/nfft;
for k=1:n_rev
    Cc(:,k)=mscohere(rev4ng(:,lr,1),rev4ng(:,lr,k),hanning(nfft),[],nfft,fs);
end
%---------------plot coherences---------
figure(70)
semilogx(fc,Cc(:,1),'b')
grid on;hold on;
for k=2:n_rev
    semilogx(fc,Cc(:,k),'r')
end
axis([0,fs/2,0,1])
xlabel('frequency [Hz]')
title('Groove Coherences')
% disp('-------------------finished--------------------') 
