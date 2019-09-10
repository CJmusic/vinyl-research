% PSD Noise analysis
% John Vanderkooy
% Jan 12, 2018
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
filename='FocusriteOpen24bit.wav';
%filename='FocusriteGen1lineXLRopen.wav';
%filename='FocusriteGen1instTSopen.wav';
%filename='DualpreOpen24bit.wav';
%filename='FocusriteOpen16bit.wav';
%filename='DualpreOpen16bit.wav';
%filename='UCA202open.wav';
%filename='UCA202short.wav';
[Focus24,fs]=audioread(filename);
Nt=length(Focus24);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-----------------choose steady portion of data----------
ns=round(5.0*fs);
nf=round(19.0*fs);
Focus24=Focus24(ns:nf,1);
%Focus16=Focus16(ns:nf,1);
%R=sqrt(12)*(rand(nf-ns+1,1)-0.5);% should plot 0dB on psd
%
R=2^-15*(rand(nf-ns+1,1)-0.5);% normal 16 bit q-noise delta^2/12
k=[1:nf-ns+1];
freq=fs/128;
%R=sin(2*pi*k*freq/fs);%test sinewave on a bin
%----------------remove integer pcm dc component----------
Idata=2^23*Focus24;%should now be all integers
dc=sum(Idata)/(nf-ns+1)
dc=round(dc);% use integer level change
Focus24=Focus24-2^-23*dc;% back to range(-1,+1)
%----------------plot histogram--------------
figure(12)
hist(Idata,512)
xlabel('q-level [24 bit]')
ylabel('number')
grid on;
%--------------get integer representation and quantize down------
nbits=22;
Focus22=2^-(nbits-1)*round(2^(nbits-1)*Focus24);
nbits=20;
Focus20=2^-(nbits-1)*round(2^(nbits-1)*Focus24);
nbits=19;
Focus19=2^-(nbits-1)*round(2^(nbits-1)*Focus24);
nbits=18;
Focus18=2^-(nbits-1)*round(2^(nbits-1)*Focus24);
nbits=17;
Focus17=2^-(nbits-1)*round(2^(nbits-1)*Focus24);
nbits=16;
Focus16=2^-(nbits-1)*round(2^(nbits-1)*Focus24);
%-----------arrange a digital noise limit plot +1/0 LSB random----
%test(:,1,1)=(round(rand(nf-ns+1,1))-0.5)/2^15; %for 16 bit system
%-----------------plot signal for each quantization-----------------------------
figure(20)
plot(Focus24(100000:200000,1),'b');
grid on;hold on;
plot(Focus22(100000:200000,1),'r');
plot(Focus20(100000:200000,1),'g');
plot(Focus19(100000:200000,1),'k');
plot(Focus18(100000:200000,1),'c');
plot(Focus17(100000:200000,1),'m');
plot(Focus16(100000:200000,1),'y');
%axis([0 1000/fs ymin ymax])
xlabel('sample #')
legend('Focus24','Focus22','Focus20','Focus19','Focus18','Focus17','Focus16');
title('low-level quantized signals')
%---------------get active portion of each test----------------
nfft=16384;
%factor=sqrt(2*fs/nfft)
[Pxx(:,1) f]=psd(Focus24(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,2) f]=psd(Focus22(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,3) f]=psd(Focus20(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,4) f]=psd(Focus19(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,5) f]=psd(Focus18(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,6) f]=psd(Focus17(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,7) f]=psd(Focus16(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,8)]=psd(R,nfft,fs,hanning(nfft),nfft/2);
%[Pxx(:,8)]=psd(R,nfft,fs);%,hanning(nfft),nfft/2);
%---------------smoothing to gauge levels----------------
for k=1:8
%Pxx(:,k)=boxsmooth(Pxx(:,k),10);
Pxx(:,k)=pwroctsmooth(Pxx(:,k),1.0);
end
%-----------------plot psd for each test-----------------------------
figure(30)
semilogx(f,10*log10(Pxx(:,1)),'b');
grid on;hold on;
semilogx(f,10*log10(Pxx(:,2)),'r');
semilogx(f,10*log10(Pxx(:,3)),'g');
semilogx(f,10*log10(Pxx(:,4)),'k');
semilogx(f,10*log10(Pxx(:,5)),'c');
semilogx(f,10*log10(Pxx(:,6)),'m');
semilogx(f,10*log10(Pxx(:,7)),'y');
semilogx(f,10*log10(Pxx(:,8)),'--k');
axis([f(2) fs/2 -110 10])
xlabel('freq[Hz]');
ylabel('PSD');
legend('Focus24','Focus22','Focus20','Focus19','Focus18','Focus17','Focus16','D2/12[16]');
title('PSD of open input preamp signals')


% disp('-------------------finished--------------------') 
