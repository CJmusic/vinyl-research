% Record Surface Noise Weighting analysis
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
filename='Silent_Spaced_1kHz.wav';
[sig,fs]=audioread(filename);
Nt=length(sig);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
figure(5)
plot(t,sig(:,1),'b') 
grid on;hold on;
plot(t,sig(:,2),'r')
axis([0 Nt/fs -.3 .3])
xlabel('time[s]')
legend('left','right')
title('All data')
%---------------------select data portion-----------------------------
ns=round(0*fs)+1;%this is where we set the analysis portion
nf=round(60*fs);
sig=sig(ns:nf,:);
lr=1;%set to 1 for left channel, 2 for right channel
Nt=length(sig);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-----------------------data plots------------
figure(10)
plot(t,sig(:,1),'b') 
grid on;hold on;
plot(t,sig(:,2),'r')
axis([0 Nt/fs -.3 .3])
xlabel('time[s]')
legend('left','right')
title('selected data')
%------test sine-----------
%sig=0.26*sin(2*pi*1000*t);
%------------------determine rms-----------
RIAArefrms=5.0;% rms cm/s, regarded as 0 dB
% if records have 40cm/s peak velocities, add 15 dB to S/N
FSref=0.262;%relative digital signal for reference peak
rms=sqrt(sum(sig.^2)/(nf-ns+1)) %still in +1-1 digital scale
dBunwtd=20*log10(rms*sqrt(2)/FSref)
%----------------------A weighting---------------
% Analog A-weighting filter according to IEC/CD 1672.
f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;
A1000 = 1.9997;
NUM = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
DEN = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
DEN = conv(conv(DEN,[1 2*pi*f3]),[1 2*pi*f2]);
%Bilinear transformation of analog design to get the digital filter. 
[b,a] = bilinear(NUM,DEN,fs);
Asig=filter(b,a,sig);
RIAArefrms=5.0;% rms cm/s, regarded as 0 dB 
FSref=0.262;%relative digital signal for reference peak
rms=sqrt(sum(Asig.^2)/(nf-ns+1)) %still in +1-1 digital scale
dBA=20*log10(rms*sqrt(2)/FSref)
%-------------------CCIR-ARM/RMS weighting----------------
Wn=2*(8400)/fs;
% good peak design Wn=8400, 2 2nd LP W=1.1Wn, 1st HP W=0.8Wn, factor=0.77*4.7544
[b,a]=butter(2,1.1*Wn);
sig1=sig;
for k=1:2
   sig1=filter(b,a,sig1);
end
%---------1xhighpass---------
[b,a]=butter(1,0.8*Wn,'high');
sig1=filter(b,a,sig1);
multfactor=0.77*4.7544;% to bring g=1 at 2kHz
CCIRsig=multfactor*sig1;
rms=sqrt(sum(CCIRsig.^2)/(nf-ns+1)) %still in +1-1 digital scale
avg=sum(abs(CCIRsig))/(nf-ns+1) %still in +1-1 digital scale
dBCCIR=20*log10(rms*sqrt(2)/FSref)
dBCCIRARM=20*log10(avg*sqrt(2)/FSref)
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
% figure(50)
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
% figure(60)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence L groove1 L groove4')
% 
% C=mscohere(rev4ng(:,1,1),rev4ng(:,1,2),hanning(nfft),[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(70)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence hann L groove1 & L groove2')
% 
% C=mscohere(rev4ng(:,1,1),rev4ng(:,2,1),hanning(nfft),[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(80)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence L groove1 & R groove1')
% 
% C=mscohere(rev4ng(:,2,1),rev4ng(:,2,2),hanning(nfft),[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(90)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence R groove1 & R groove2')

% disp('-------------------finished--------------------') 
