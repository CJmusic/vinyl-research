%CCIR filter reverse engineered
clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
disp('CCIR/ARM filter designs')

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
fs=44100;
time=([1:N]'-1)/fs; % make time positive column vector starting at zero
signal=zeros(N,1);
impulse=zeros(N,1);
impulse(1)=1;%unit impulse
%--------
Wn=2*(8400)/fs;
% fair design Wn=8700, 2nd LP, 1st HP factor 1.05*4.7544
% good 44.1 peak design Wn=8400, 2 2nd LP W=1.1Wn, 2 (1st LP W=1.2Wn ?), 1st HP W=0.8Wn, factor=0.77*4.7544
%-------3 2nd-order LP + one 1st-order LP-------------
[b,a]=butter(2,1.3*Wn);
signal=impulse;
for k=1:1
   signal=filter(b,a,signal);
end
[b,a]=butter(1,1.4*Wn);
for k=1:2
   signal=filter(b,a,signal);
end
%---------1st order highpass---------
[b,a]=butter(1,0.8*Wn,'high');
signal=filter(b,a,signal);
multfactor=0.77*4.7544;% to bring gain=1 at 2kHz
signal=multfactor*signal;
%--------------------time plot----------------------------
figure(10)
plot(time,signal,'b');
grid on;
axis([0 .001 -1 1]);
title('CCIR/ARM 44.1kHz rev eng impulse signal')
xlabel('Time [s]')
%-----------------freq plot-----------------
OUTPUT=fft(signal);
f=([1:N/2+1]'-1)*fs/N;
figure(20);
semilogx(f,20*log10(abs(OUTPUT(1:floor(N/2+1)))),'b');
grid on;hold on;
semilogx(fr,CCIR,'r');
pk=10*ceil(max(20*log10(abs(OUTPUT)))/10);
axis([5*fs/N,fs/2,pk-60,pk])
legend('CCIR rev engineered','CCIR/ARM table','Location','Best');%,'ITUMAG','ITUFILT','ITUFILT2','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('magnitude response');
%--------------phase response--------------------
% figure(30);
% semilogx(f,(180/pi)*angle(OUTPUT(1:floor(N/2+1))),'b');
% grid on;
% axis([fs/N,fs/2,-180,180])
% xlabel('Frequency [Hz]')
% ylabel('phase')
% title('net phase response');
%------------------yule-walker 44.1kHz filter design----------------
fs=44100;
frdc=[0 31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000  fs/2];
CCIR=[-100 -35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8 -32];
Wn=2*frdc/fs;
CCIRmag=10.^(CCIR/20);
[b,a]=yulewalk(12,Wn,CCIRmag);
[d,c]=butter(1,2*370/fs,'high');% this corrects DC-LF with highpass
%------------------compactify into single filter----------------
%signal=filter(conv(b,d),conv(a,c),impulse);
fb=conv(b,d);ea=conv(a,c);
signal=filter(fb,ea,impulse);
%-----------------plot time---------------
% time=([1:N]'-1)/fs;
% figure(40)
% plot(time,signal,'b');
% grid on;
% axis([0 .001 -1 1]);
% title('CCIR/ARM julewalk 44.1k impulse signal')
% xlabel('Time [s]')
%----------------------plot CCIR yulewalk response------------
f=[0:N/2]*fs/N;
OUTPUT=fft(signal);
%--------------------------------------------------------------------
figure(50);
semilogx(f,20*log10(abs(OUTPUT(1:floor(N/2+1)))),'b');
grid on;hold on;
semilogx(frdc,CCIR,'r');
pk=10*ceil(max(20*log10(abs(OUTPUT)))/10);
axis([fs/N,fs/2,pk-60,pk])
legend('CCIR yulewalk','CCIR/ARM table','Location','Best');%,'ITUMAG','ITUFILT','ITUFILT2','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('44.1 CCIR magnitude response');
%------------------yule-walker 96kHz filter design----------------
fs=96000;
frdc=[0 31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000 25000 30000 fs/2];
CCIR=[-inf -35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8 -35 -50 -inf];
Wn=2*frdc/fs;
CCIRmag=10.^(CCIR/20);
[b,a]=yulewalk(12,Wn,CCIRmag);%%%%%%%%%%%%%%%%%%%%%
[d,c]=butter(1,2*750/fs,'high');% this corrects DC-LF with highpass
%------------------compactify into single filter----------------
%signal=filter(conv(b,d),conv(a,c),impulse);
fb=conv(b,d);ea=conv(a,c);
signal=filter(fb,ea,impulse);
%-----------------plot time---------------
% time=([1:N]'-1)/fs;
% figure(60)
% plot(time,signal,'b');
% grid on;
% axis([0 .001 -1 1]);
% title('CCIR/ARM julewalk 96k impulse signal')
% xlabel('Time [s]')
%----------------------plot yulewalk freq------------
f=[0:N/2]*fs/N;
OUTPUT=fft(signal);
%--------------------------------------------------------------------
figure(70);
semilogx(f,20*log10(abs(OUTPUT(1:floor(N/2+1)))),'b');
grid on;hold on;
semilogx(frdc,CCIR,'r');
pk=10*ceil(max(20*log10(abs(OUTPUT)))/10);
axis([fs/N,fs/2,pk-60,pk])
legend('CCIR yulewalk','CCIR/ARM table','Location','Best');%,'ITUMAG','ITUFILT','ITUFILT2','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('96k CCIR magnitude response');
%----------------------------------------------------------


disp('-----------------------------finished------------------------')


