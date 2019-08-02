% Record Surface Noise analysis
% using xcorrelation between groove signals
% John Vanderkooy
% Jan. 2019
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
filename='/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/vinyl-research/matlab_code/ChirpSilentLeadin.wav'; %filename='leadin35sec.wav';
[sig,fs]=audioread(filename);
ns=round(35.0*fs);%this is where we set the analysis portion
nf=round(71.0*fs);
sig=sig(ns:nf,:);
Nt=length(sig);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-----------------------data plots------------
figure(12)
plot(t,sig(:,1),'b') 
grid on;hold on;
plot(t,sig(:,2),'r')
axis([0 Nt/fs -.1 .1])
xlabel('time[s]')
legend('left','right')
title('All data')
%------------------generate separate groove results-----------
% we'll make n_rev files
n_per_rev=round(fs*60/33.3333)
n_rev=floor(Nt/n_per_rev) % number of grooves

% choose the data to plot
for ng=1:n_rev;% n_rev % groove number
    rev4ng(:,:,ng)=sig(1+(ng-1)*n_per_rev:ng*n_per_rev,:);
end
%---------------correlation between grooves--------------
for ng=1:n_rev
    cf(ng)=sum(rev4ng(:,1,1).*rev4ng(:,1,ng),1);
end
figure(13)
plot(cf)
xlabel('groove number')
ylabel('correlation coefficient')
grid on
%---------time variable in plots------------
tg=t(1:n_per_rev);
minplot=-0.1;
maxplot=+0.1;
%-----------------plot signal for each groove-----------------------------
figure(14)
plot(tg,rev4ng(:,1,1),'b') 
grid on;hold on;
plot(tg,rev4ng(:,1,2),'r'); % variables are sample#, left/right, groove#
plot(tg,rev4ng(:,1,3),'g') 
plot(tg,rev4ng(:,1,4),'k') 
plot(tg,rev4ng(:,1,5),'y');
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
% figure(18)
% plot(tg,rev4ng(:,1,1),'b') 
% grid on;hold on;
% plot(tg,rev4ng(:,2,1),'r') 
% axis([0 n_per_rev/fs minplot maxplot])
% xlabel('time[s]')
% legend('groove1 left','groove1 right')
% title('Declicked stereo groove')
% 
% %---------------cross correlation properties-----------
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
% %----------------------coherence-----------------
% nfft=1024;
% fc=([0:nfft/2])*fs/nfft;
% 
% % rev4ng(:,1,1)=rand(n_per_rev,1);%testing coherence
% % rev4ng(:,1,4)=rand(n_per_rev,1);
% 
% C=mscohere(rev4ng(:,1,1),rev4ng(:,1,4),hanning(nfft),[],nfft,fs);
% %C=mscohere(rev4ng(:,1,1),rev4ng(:,2,1),[],[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(25)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence L groove1 L groove4')
% 
% C=mscohere(rev4ng(:,1,1),rev4ng(:,1,2),hanning(nfft),[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(26)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence hann L groove1 & L groove2')
% 
% C=mscohere(rev4ng(:,1,1),rev4ng(:,2,1),hanning(nfft),[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(27)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence L groove1 & R groove1')
% 
% C=mscohere(rev4ng(:,2,1),rev4ng(:,2,2),hanning(nfft),[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(28)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence R groove1 & R groove2')

% disp('-------------------finished--------------------') 
