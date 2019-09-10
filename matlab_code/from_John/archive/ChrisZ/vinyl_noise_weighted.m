% Record Surface Noise Weighting analysis
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
%---------------calibration constants-------------------
%the following reference level must be measured from a trusted test disc.
FSref=0.262;%measured relative digital signal for peak (5cm/s rms=7cm/s peak)
RIAArefrms=5.0;% rms cm/s, regarded as 0 dB.  We want signal in these units
Vpeak=40;% 40 cm/s is regarded as peak signal
% peak_signal/ref_signal=Vpeak/(RIAArefrms*sqrt(2))-> 15.05 dB
%----------------calculate signal amplitude factor---------
sigfactor=sqrt(2)/FSref;% sig normalize factor so that unity is 5rms/7peak 
%---------------------------------------------
lr=2;%left=1, right=2
filename='J15_tr5_outer.wav';ts=0;tf=21.5;
filename='J15_tr5_inner.wav';ts=0;tf=21.5;
%filename='DualpreOpen16bit.wav';ts=0;tf=20;
%filename='Silent_Spaced_1kHz.wav';ts=0;tf=60;
%filename='Silent_Spaced_1kHz.wav';ts=74;tf=97;
%filename='Silent_Spaced_1kHz.wav';ts=103;tf=118;% reference 5cm/s rms
%filename='RCArequiem1.wav';ts=4;tf=8;
%filename='RCArequiem2.wav';ts=4;tf=14;
%filename='RCArequiem3.wav';ts=4;tf=14;
%filename='RCArequiem4.wav';ts=6;tf=14;
%filename='SMCstart1.wav';ts=2.5;tf=9;
%filename='SMCstart2.wav';ts=2.2;tf=6.5;
%filename='Joni1.wav';ts=3.5;tf=6.0;
%filename='Joni2.wav';ts=2;tf=5.3;
%filename='DGG1.wav';ts=3;tf=9;
%filename='DGG2.wav';ts=2.5;tf=9.5;
%filename='analogmagik_silent_1kHz.wav';ts=0;tf=8.5;
[sig,fs]=audioread(filename);
Nt=length(sig);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
figure(5)
plot(t,sig(:,1),'b') 
grid on;hold on;
plot(t,sig(:,2),'r')
axis([0 Nt/fs ylim])
xlabel('time[s]')
legend('left','right')
title('All data')
%---------------------select data portion-----------------------------
ns=round(ts*fs)+1;%this is where we set the analysis portion
nf=round(tf*fs);
sig=sigfactor*sig(ns:nf,:);% now normalized in units of RIAArefrms*sqrt(2)
Nt=length(sig);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=(0:Nt-1)/fs;%column vector
%-----------------------data plots------------
figure(10)
plot(t,sig(:,1),'b') 
grid on;hold on;
plot(t,sig(:,2),'r')
axis([0 t(Nt) ylim])
xlabel('time[s]')
ylabel('amplitude')
legend('left','right')
title('selected unweighted data')
%------test sine and impulse----------
%sig=zeros(Nt,1);sig(1000)=2^14;ns=1;nf=Nt;
%sig=FSref*sin(2*pi*1000*t);% should be 0 dB
%------------------determine rms unweighted-----------
rms=sqrt(sum(sig.^2)/(nf-ns+1)); %should be unity for ref signal
dBunwtd=20*log10(rms);
disp(['rms:' num2str(rms) ' dBunwtd:' num2str(dBunwtd)])
%-----------------filter to enhance clicks?-------------
[b,a]=butter(1,0.7,'high');
sigf=filtfilt(b,a,sig);
figure(15)
plot(t,sigf(:,1),'b') 
grid on;hold on;
plot(t,sigf(:,2),'r')
axis([xlim ylim])
xlabel('time[s]')
ylabel('amplitude')
legend('left','right')
title('selected filtered data')
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
sigA=filter(b,a,sig);
rms=sqrt(sum(sigA.^2)/(nf-ns+1)); %still in +1-1 digital scale
dBA=20*log10(rms);
disp(['rmsA:' num2str(rms) ' dBA:      ' num2str(dBA)])
%-------------------CCIR-ARM/RMS weighting----------------
Wn=2*(8400)/fs;
% good peak design Wn=8400, 2 2nd LP W=1.1Wn, 1st HP W=0.8Wn, factor=0.77*4.7544
[b,a]=butter(2,1.1*Wn);
sigC=sig;
for k=1:2
   sigC=filter(b,a,sigC);%recursive use
end
%---------1xhighpass---------
[b,a]=butter(1,0.8*Wn,'high');
sigC=filter(b,a,sigC);
multfactor=0.77*4.7544;% to bring gain=1 at 2kHz
sigC=multfactor*sigC;
rms=sqrt(sum(sigC.^2)/(nf-ns+1)); %still in +1-1 digital scale
dBCCIRrms=20*log10(rms);
avg=sum(abs(sigC))/(nf-ns+1); %avg responding, not sine corrected (-1.1dB)
dBCCIRarm=20*log10(avg);
disp(['rmsC:' num2str(rms) ' dBCCIRrms:' num2str(dBCCIRrms)])
disp(['avgC:' num2str(avg) ' dBCCIRavg:' num2str(dBCCIRarm)])
%--------------------plot filtered signals------------
figure(20)
plot(t,sigA(:,lr),'b') 
grid on;hold on;
plot(t,sigC(:,lr),'r')
%axis([: : -.3 .3])
xlabel('Time[s]')
legend('Awtd','CCIR')
title('selected weighted data')
%---------------------spectrum plots--------------------
nfft=16384;
%factor=sqrt(2*fs/nfft)
[Pxx(:,lr),f]=pwelch(sig(:,lr),hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
[Pxx2(:,lr),f]=pwelch(sigA(:,lr),hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
[Pxx3(:,lr),f]=pwelch(sigC(:,lr),hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
%---------------smoothing if necessary----------------
% Pxx(:,1)=pwroctsmooth(Pxx(:,1),0.1);% these need to have full complex Pxx
% Pxx2(:,1)=pwroctsmooth(Pxx2(:,1),0.1);
% Pxx3(:,1)=pwroctsmooth(Pxx3(:,1),0.1);
%-----------------plot psd for each test-----------------------------
figure(30)
semilogx(f,10*log10(Pxx(:,lr)),'b');
grid on;hold on
semilogx(f,10*log10(Pxx2(:,lr)),'r');
semilogx(f,10*log10(Pxx3(:,lr)),'g');
axis([f(2) fs/2 -140 -40])
xlabel('freq[Hz]');
ylabel('Power Spectrum');
legend('raw psd' ,'A-wtd','CCIR wtd');
title('PSD')
disp('-------------------finished--------------------') 
