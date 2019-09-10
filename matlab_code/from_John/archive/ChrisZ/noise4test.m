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
filename='FocusriteOpen24bit.wav';
[Focus24,fs]=audioread(filename); 
filename='FocusriteOpen16bit.wav';
[Focus16,fs]=audioread(filename); 
filename='DualpreOpen24bit.wav';
[Dualpre24,fs]=audioread(filename); 
filename='DualpreOpen16bit.wav';
[Dualpre16,fs]=audioread(filename); 
Nt=length(Focus24);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-----------------choose steady portion of data----------
ns=round(1.0*fs);
nf=round(19.0*fs);
Focus24=Focus24(ns:nf,1);
Focus16=Focus16(ns:nf,1);
Dualpre24=Dualpre24(ns:nf,1);
Dualpre16=Dualpre16(ns:nf,1);
%-----------arrange a digital noise limit plot +/0 LSB random----
test(:,1,1)=(round(rand(nf-ns+1,1))-0.5)/2^15; %for 16 bit system

ymin=-1e-2;
ymax=+1e-2;
%-----------------plot signal for each test-----------------------------
figure(20)
plot(Focus24(10000:100000,1),'b');
grid on;hold on;
plot(Focus16(10000:100000,1),'r');
plot(Dualpre24(10000:100000,1),'g');
plot(Dualpre16(10000:100000,1),'k'); 
plot(test(1:1000,1),'c'); 
%axis([0 1000/fs ymin ymax])
xlabel('time[s]')
legend('Focus24','Focus16','Dualpre24','Dualpre16')
title('low-level signals')
%---------------get active portion of each test----------------
nfft=16384;
[Pxx(:,1) f]=psd(Focus24(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,2) f]=psd(Focus16(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,3) f]=psd(Dualpre24(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,4) f]=psd(Dualpre16(:,1),nfft,fs,hanning(nfft),nfft/2);
[Pxx(:,5) f]=psd(test,nfft,fs,hanning(nfft),nfft/2);
%-----------------plot psd for each test-----------------------------
figure(30)
semilogx(f,10*log10(Pxx(:,1)),'b');
grid on;hold on;
semilogx(f,10*log10(Pxx(:,2)),'r');
semilogx(f,10*log10(Pxx(:,3)),'g');
semilogx(f,10*log10(Pxx(:,4)),'k');
semilogx(f,10*log10(Pxx(:,5)),'c');
axis([f(2) fs/2 ylim])
xlabel('freq[Hz]');
ylabel('PSD');
legend('Focus24','Focus16','Dualpre24','Dualpre16','1 LSB 16');
title('PSD of open input preamp signals')


% disp('-------------------finished--------------------') 
