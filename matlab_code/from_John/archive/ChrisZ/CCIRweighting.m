%CCIR filter reverse engineered
clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

Po=2e-5;% SPL ref
try
    pkg load signal; %for Octave
catch
end
% -----------------------------------------------------------------
N=2^14;
%----------------------CCIR/ARM dB table-------------------------------
fr=[31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000];
CCIR=[-35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8];
%------------------this 44.1kHz fit needs optimization---------------------
fs=44100;disp(['Total duration [s] ' num2str(N/fs)])
time=([1:N]'-1)/fs; % make time positive column vector starting at zero
signal=zeros(N,1);
impulse=zeros(N,1);
impulse(1)=1;%unit impulse

figure(10)
plot(time,impulse,'b');
grid on;
title('input signal')
xlabel('Time [s]')
Wn=2*(8400)/fs;
% fair design Wn=8700, 2nd LP, 1st HP factor 1.05*4.7544
% good 44.1 peak design Wn=8400, 2 2nd LP W=1.1Wn, 2 (1st LP W=1.2Wn ?), 1st HP W=0.8Wn, factor=0.77*4.7544
%-------3 repeated 2nd-order LP + one 1st-order LP-------------
[b,a]=butter(2,1.2*Wn);
signal=impulse;
for k=1:2
   signal=filter(b,a,signal);
end
[b,a]=butter(2,1.2*Wn);
for k=1:1
   signal=filter(b,a,signal);
end
%---------1xhighpass---------
[b,a]=butter(1,1.12*Wn,'high');
signal=filter(b,a,signal);
multfactor=1.0*4.7544;% to bring g=1 at 2kHz
signal=multfactor*signal;
% %------------------this is 96kHz optimization---------------------
fs=96000;disp(['Total duration [s] ' num2str(N/fs)])
% time=([1:N]'-1)/fs; % make time positive column vector starting at zero
% signal=zeros(N,1);
% impulse=zeros(N,1);
% impulse(1)=1;%unit impulse
% 
% figure(20)
% plot(time,impulse,'b');
% grid on;
% title('input signal')
% xlabel('Time [s]')

% Wn=2*(10700)/fs;
% %-------3 repeated 2nd-order LP + one 1st-order LP-------------
% [b,a]=butter(2,Wn);
% signal=impulse;
% for k=1:3
%    signal=filter(b,a,signal);
% end
% [b,a]=butter(1,Wn);
% for k=1:1
%    signal=filter(b,a,signal);
% end
% %---------1xhighpass---------
% [b,a]=butter(1,Wn,'high');
% signal=filter(b,a,signal);
% multfactor=5.7544;% to bring g=1 at 2kHz
% signal=multfactor*signal;
%---------------------------
OUTPUT=fft(signal);
%-------------------------plots----------------------------
figure(30)
plot(time,signal,'b');
grid on;
title('CCIR/ARM impulse signal')
xlabel('Time [s]')

figure(40)
plot(time,signal,'b');
grid on;
axis([0 .001 -1 1]);
title('CCIR/ARM impulse signal')
xlabel('Time [s]')
%---------------frequency domain----------------------
%fftfactor=sqrt(2/(N*fs)); %Pa/root_Hz single sided
f=([1:N/2+1]'-1)*fs/N;


%--------------------------------------------------------------------
figure(50);
semilogx(f,20*log10(abs(OUTPUT(1:floor(N/2+1)))),'b');
grid on;hold on;
semilogx(fr,CCIR,'r');
pk=10*ceil(max(20*log10(abs(OUTPUT)))/10);
axis([35*fs/N,0.3*fs/2,pk-20,pk])
legend('CCIR rev engineered','CCIR/ARM table','Location','Best');%,'ITUMAG','ITUFILT','ITUFILT2','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('magnitude response');

% figure(55);
% semilogx(f,(180/pi)*angle(OUTPUT(1:floor(N/2+1))),'b');
% grid on;
% axis([fs/N,fs/2,-180,180])
% xlabel('Frequency [Hz]')
% ylabel('phase')
% title('net phase response');

%------------------yule-walker filter design----------------
frdc=[0 31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000 fs/2];
CCIR=[-300 -35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8 -300];
Wn=2*frdc/fs;
CCIRmag=10.^(CCIR/20);
[b,a]=yulewalk(8,Wn,CCIRmag);
signal=impulse;
signal=filter(b,a,signal);
[d,c]=butter(1,2*1070/fs,'high');% this corrects DC with LF highpass
signal=filter(d,c,signal);
%----------------------plot yulewalk & CCIR table------------
%fftfactor=sqrt(2/(N*fs)); %Pa/root_Hz single sided
f=([1:N/2+1]'-1)*fs/N;
OUTPUT=fft(signal);
%--------------------------------------------------------------------
figure(60);
semilogx(f,20*log10(abs(OUTPUT(1:floor(N/2+1)))),'b');
grid on;hold on;
semilogx(frdc,CCIR,'r');
pk=10*ceil(max(20*log10(abs(OUTPUT)))/10);
axis([fs/N,fs/2,pk-60,pk])
legend('CCIR yulewalk','CCIR/ARM table','Location','Best');%,'ITUMAG','ITUFILT','ITUFILT2','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('magnitude response');

% %--------------------ITU-R 468 parameterization--------------------
% b0=1.2463e-4;
% a6=-4.7373e-24;a5=1.3066e-19;a4=2.0438e-15;a3=-2.1182e-11;a2=-1.3639e-7;a1=5.5595e-4;
% 
% f2=2000;% for CCIR/ARM
% h1=a6*f2^6+a4*f2^4+a2*f2^2+1;
% h2=a5*f2^5+a3*f2^3+a1*f2;
% ITUMAG2=b0*f2/sqrt(h1^2+h2^2);%value of ITU filter gain at 2kHz
% 
% h1=a6*f.^6+a4*f.^4+a2*f.^2+1;
% h2=a5*f.^5+a3*f.^3+a1*f;
% ITUMAG=(1/ITUMAG2)*b0*f./sqrt(h1.^2+h2.^2);
% %----------transform to Matlab form---problem!!!-------
% b0=1.2463e-4;%s-plane numerator
% a=[-a6 a5 a4 -a3 -a2 a1 1];%s-plane denominator, now all positive
% b=[b0 0 0 0 0 0 0];
% [b,a]=bilinear(b,a,fs)
% itufilt=1e-24*filter(b,a,impulse);
% figure(55);
% plot(itufilt)
% 
% ITUFILT=fft(itufilt);
% figure(70)
% plot(2*pi*f,20*log10(abs(ITUFILT(1:N/2+1))))
% %this s-plane polynomial seems to use omega=2Pi f, so plot below is
% %modified.  There seem to be small errors, too high for roundoff?
% 
% %----------transform to s-plane first to avoid roundoff?!-------
% b0=1.2463e-4;%s-plane numerator
% pi1=2*pi;pi2=2*pi*pi1;pi3=2*pi*pi2;pi4=2*pi*pi3;pi5=2*pi*pi4;pi6=2*pi*pi5;
% a=[-a6/pi6 a5/pi5 a4/pi4 -a3/pi3 -a2/pi2 a1/pi1 1];%s-plane denominator, now all+
% b=[b0 0 0 0 0 0 0];
% [b,a]=bilinear(b,a,fs)
% itufilt2=pi6*1e-24*filter(b,a,impulse);
% figure(75);
% plot(itufilt2)
% grid on
% 
% ITUFILT2=fft(itufilt2);
% figure(80)
% semilogx(f,-400+20*log10(abs(ITUFILT2(1:N/2+1))))
% grid on
% %this s-plane polynomial seems to use omega=2Pi f, so plot is modified.
% %There seem to be small errors, too high for roundoff?


disp('-----------------------------finished------------------------')


