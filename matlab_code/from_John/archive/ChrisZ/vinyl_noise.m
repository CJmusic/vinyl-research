% Record Surface Noise analysis
% John Vanderkooy
% March 2018
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
filename='Normalize1k2_4V.wav';
[N2_4,fs]=audioread(filename); 
Nt=length(N2_4)
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%------------------------------------------------
figure(10)
plot(t,N2_4(:,1),'b') 
grid on;hold on;
plot(t,N2_4(:,2),'r')
axis([0 Nt/fs -1 1])
xlabel('time[s]')
legend('left','right')
title('Original')
%-----------------------system noise------------
filename='Silence_no_record.wav';
[silnr,fs]=audioread(filename); 
Nt=length(silnr)
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%------------------------------------------------
figure(12)
plot(t,silnr(:,1),'b') 
grid on;hold on;
plot(t,silnr(:,2),'r')
axis([0 Nt/fs -1 1])
xlabel('time[s]')
legend('left','right')
title('Original')

%-----------------------silent track-----------
filename='Dec20_silent_track.wav';
[sil,fs]=audioread(filename); 
Nt=length(sil)
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%------------------------------------------------
figure(14)
plot(t,sil(:,1),'b') 
grid on;hold on;
plot(t,sil(:,2),'r')
axis([0 Nt/fs -1 1])
xlabel('time[s]')
legend('left','right')
title('Original')

%-------------------------------
filename='Dec20_40cmps.wav';
[D40cs,fs]=audioread(filename); 
Nt=length(D40cs)
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%------------------------------------------------
figure(16)
plot(t,D40cs(:,1),'b') 
grid on;hold on;
plot(t,D40cs(:,2),'r')
axis([0 Nt/fs -1 1])
xlabel('time[s]')
legend('left','right')
title('Original')

%---------------------select analysis segment-------------
starttime=1.0;
ns=round(starttime*fs)+1;
N=2^16;%length of analysis
disp(' ')
disp(['N: ' num2str(N) ' duration: ' num2str(N/fs)])
t=linspace(0,(N-1)/fs,N)';%column vector
sil=sil(ns:ns+N-1,:);
silnr=silnr(ns:ns+N-1,:);

%----------------------------------------
figure(20)
plot(t,sil(:,1),'b') 
grid on;hold on;
plot(t,sil(:,2),'r')
axis([0 N/fs -.1 .1])
xlabel('time[s]')
legend('silent track','system noise')
title('silent track')
%----------------------------------------
figure(25)
plot(t,silnr(:,1),'b') 
grid on;hold on;
plot(t,silnr(:,2),'r')
axis([0 N/fs -.01 .01])
xlabel('time[s]')
legend('silent track','system noise')
title('noise:no record')

%-------------------spectral plots-----------------------
f=([1:N/2+1]-1)*fs/N;
f(1)=1e-10;% prevents log plot warning in Octave
SIL=fft(sil(:,2));
SILNR=fft(silnr(:,2));

figure(100);
semilogx(f,20*log10(abs(SIL(1:N/2+1))),'b');
grid on;hold on;
semilogx(f,20*log10(abs(SILNR(1:N/2+1))),'r');
axis([1 fs/2 -60 40])
legend('silent track','system noise')
xlabel('frequency [Hz]');ylabel('dB')
title(['Spectra']);
%---------------------------------------
SILm=pwroctsmooth(SIL,0.1);
SILNRm=pwroctsmooth(SILNR,0.1);

figure(200);
semilogx(f,20*log10(abs(SILm(1:N/2+1))),'b');
grid on;hold on
semilogx(f,20*log10(abs(SILNRm(1:N/2+1))),'r');
axis([1 fs/2 -60 40])
legend('silent track','system noise')
xlabel('frequency [Hz]');ylabel('dB')
title(['Smoothed Spectra']);

disp('-------------------finished--------------------') 
