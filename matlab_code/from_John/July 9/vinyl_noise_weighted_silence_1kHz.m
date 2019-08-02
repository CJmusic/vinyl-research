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
% reference level must be measured from the test disc 1kHz track Vp=7.0cm/s.
% 40 cm/s sine is regarded as peak signal, 0dB (28.3 cm/s rms)
% peak_signal/ref_signal=40/7 -> 15.05 dB
%---------------------------------------------
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
% path_folder = strcat('/Volumes/AUDIOBANK/audio_files/pressings/', folder, '/')
% filename = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/CZr9pioneer.wav');ts=0;tf=97;tts=108;ttf=135;
filename = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r26-96khz-declicked.wav');ts=0;tf=97;tts=108;ttf=135;
lr=1;%left=1, right=2
%filename='r9technics_declick_sil&1kHz.wav';ts=0;tf=92;tts=105;ttf=145;
% filename='r9pioneer_declick_sil&1kHz.wav';ts=0;tf=97;tts=103;ttf=120;
%filename='28b_declick_full_transition.wav';ts=10.0;tf=106.5;tts=115;ttf=129.0;
%filename='052319_pioneer3150.wav';ts=0;tf=20;tts=30;ttf=70;
%filename='lacquerA_sil_1kHz_beginning.wav';ts=0;tf=13;tts=18;ttf=24;
%filename='lacquerA_silence_edited_onto_1kHz_declicked.wav';ts=0;tf=40;tts=50;ttf=80;
%filename='lacquerA_silence_edited_onto_1kHz.wav';ts=0;tf=40;tts=50;ttf=80;
%filename='lacquerA_silence_1kHz.wav'; ts=0;tf=40;tts=50;ttf=90;
%filename='28b_silence_1kHz.wav';ts=0;tf=94;tts=105;ttf=145;
%filename='analogmagik_silent_1kHz_96k.wav';ts=0;tf=8.5;tts=8.7;ttf=15.7;%96
[sig,fs]=audioread(filename);
Nt=length(sig);
disp([num2str(filename) ' fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%----------------------
figure(5)
plot(t,sig(:,1),'b') 
grid on;hold on;
plot(t,sig(:,2),'r')
axis([0 Nt/fs ylim])
xlabel('time[s]')
legend('left','right')
title('All raw data')
%-------------------select reference portion-----------------
nts=round(tts*fs);
ntf=round(ttf*fs);
ref=sig(nts:ntf,1);% this is 7 cm/sec reference, left channel


rmsref=sqrt(sum(ref.^2)/(ntf-nts+1)); % for 7 cm/s peak
%% RMS VALUE weighted bynumber of points - CZ 



disp(['rmsref raw: ' num2str(rmsref)])


%------------signal amplitude factor for 0dB=40cm/s peak---------

peakref=sqrt(2)*rmsref*40/7; %digital value of peak level





%---------------------select data portion-----------------------------
ns=round(ts*fs)+1;%this is where we set the analysis portion
nf=round(tf*fs);


sig=sig(ns:nf,:)/peakref;% now normalized to 40cm/s peak


Nt=length(sig);
disp(['noise  fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=(0:Nt-1)'/fs;%column vector
%------------prefilter signal------

%% REMOVE DC w/ High Pass at 1.5 Hz
[b,a]=butter(2,2*1.5/fs,'high');%CD HP. Choosing 20Hz removes arm resonance

sig=filter(b,a,sig);
% sig(:,1)=rand(Nt,1)-0.5;
% sig(:,1)=sqrt(2)*sin(2*pi*1000*t);%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%
% sig(:,2)=sqrt(2)*sin(2*pi*2000*t);
%-----------------------data plots------------
figure(10)
plot(t,sig(:,1),'b') 
grid on;hold on;
plot(t,sig(:,2),'r')
axis([0 t(Nt) ylim])
xlabel('time[s]')
ylabel('amplitude')
legend('left','right')
title('selected data')
%------------------determine rms unweighted-----------
rms=sqrt(sum(sig.^2)/(nf-ns+1)); %this should be -15dB for ref signal
dBunwtd=20*log10(rms);
disp(['unweighted sig rms:' num2str(rms)  '  dBunwtd:' num2str(dBunwtd)])
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
title('click enhanced filtered data')
%----------------------A weighting 96kHz---------------
% Analog A-weighting filter according to IEC/CD 1672.
f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
%f4 = 12194.217;%for infinite fs
f4 = 14100;%for fs=96kHz, makes magnitude -9.5dB at 20kHz
A1000 = 1.9997;
NUM = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
DEN = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
DEN = conv(conv(DEN,[1 2*pi*f3]),[1 2*pi*f2]);
%Bilinear transformation of analog design to get the digital filter. 
[b,a] = bilinear(NUM,DEN,fs);
sigA=filter(b,a,sig);
%------------prepare fractional detector----------
fracexp=1.5; disp(['fractional exponent: ' num2str(fracexp)])
x=pi*[0:1000]/1000;
s=(sum(sin(x).^fracexp)/1000)^(1/fracexp);%integrate over 1/2 period
factor=1/(sqrt(2)*s);% this gives 1/sqrt2 for a unit amplitude sinewave 
%------------------------------------
rmsA=sqrt(sum(sigA.^2)/(nf-ns+1)); 
dBA=20*log10(rmsA);
disp(['rmsA:' num2str(rmsA) '   dBArms:   ' num2str(dBA)])
fracA=(sum(abs(sigA).^fracexp)/(nf-ns+1)).^(1/fracexp);
dBAfrac=20*log10(fracA*factor);
disp(['fracA: ' num2str(fracA) '  dBAfrac: ' num2str(dBAfrac)])
avgA=sum(abs(sigA)*pi/(2*sqrt(2)))/(nf-ns+1);
dBAavg=20*log10(avgA);
disp(['avgA:' num2str(avgA) '   dBAavg:   ' num2str(dBAavg)])
%------------------yule-walker 96kHz CCIRarm filter design----------------
if fs ~= 96000;
   disp('incompatible sampling freq')
   halt
end
frdc=[0 31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000 25000 30000 fs/2];
CCIR=[-inf -35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8 -35 -50 -inf];
Wn=2*frdc/fs;
CCIRmag=10.^(CCIR/20);
[b,a]=yulewalk(12,Wn,CCIRmag);%%%%%%%%%%%%%%%%%%%%%
[d,c]=butter(1,2*750/fs,'high');% this corrects DC-LF with highpass
%----compactify into single filter---
fb=conv(b,d);ea=conv(a,c);
sigC=filter(fb,ea,sig);
%--------------------dB calcs--------------------------
rms=sqrt(sum(sigC.^2)/(nf-ns+1)); 
dBCCIRrms=20*log10(rms);
disp(['rmsC:' num2str(rms) '   dBCCIRrms:' num2str(dBCCIRrms)])
fracC=(sum(abs(sigC).^fracexp)/(nf-ns+1)).^(1/fracexp);
dBCfrac=20*log10(fracC*factor);
disp(['fracC: ' num2str(fracC) '  dBCIRfrac: ' num2str(dBCfrac)])
avgC=sum(abs(sigC)*pi/(2*sqrt(2)))/(nf-ns+1); %avg responding, sine corrected
dBCCIRarm=20*log10(avgC);
disp(['avgC:' num2str(avgC) '   dBCIRavg:' num2str(dBCCIRarm)])
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
nfft=2^14;
%factor=sqrt(2*fs/nfft)
[Pxx(:,lr),f]=pwelch(sig(:,lr),hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
[Pxx2(:,lr),f]=pwelch(sigA(:,lr),hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
[Pxx3(:,lr),f]=pwelch(sigC(:,lr),hanning(nfft,'periodic'),nfft/2,nfft,fs,'psd');
%-----------------plot psd for each test-----------------------------
figure(30)
semilogx(f,10*log10(Pxx(:,lr)),'b');
grid on;hold on
semilogx(f,10*log10(Pxx2(:,lr)),'r');
semilogx(f,10*log10(Pxx3(:,lr)),'g');
axis([f(2) fs/2 -160 -60])
xlabel('freq[Hz]');
ylabel('Power Spectrum');
legend('raw psd' ,'A-wtd','CCIR wtd');
title([filename ' PSD'])
disp('-------------------finished--------------------') 
