%A-weighting digital filter
clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
try
    pkg load signal; %for Octave
catch
end

addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')


%--------------------------------------
fs=96000; 
N=2^14;disp(['Total duration [s] ' num2str(N/fs)])
time=([1:N]'-1)/fs; % make time positive column vector starting at zero
signal=zeros(N,1);
signal(10)=1;%unit impulse
%-------------------------------------------------------
output=audio_Aweighting(signal);

figure(30)
plot(time,output,'b');
grid on;
axis([0 .003 -.1 .1]);
title('A-wtd impulse signal')
xlabel('Time [s]')
%---------------frequency domain----------------------
f=([1:N/2+1]'-1)*fs/N;
OUTPUT=fft(output);

%% -------------------------A-wtg table-----------------------------------
fr=[0 6.3 10 16 25 40 63 100 160 250 400 630 1000 1600 2500 4000 6300 10000 16000 20000 22050];
dBAtable=[-inf -85.4 -70.4 -56.7 -44.7 -34.6 -26.2 -19.1 -13.4 -8.6 -4.8 -1.9 0 1.0 1.3 1.0 -0.1 -2.5 -6.6 -9.3 -11.];
%---------------------------------------

figure(60);
semilogx(fr,dBAtable,'k');
grid on;hold on
semilogx(f,20*log10(abs(OUTPUT(1:floor(N/2+1)))),'b');
pk=10*ceil(max(20*log10(abs(OUTPUT)))/10);
axis([fs/N,fs/2,pk-60,pk])
legend(' A-wtg table','A-wtg 96 response','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('bilinear magnitude response');
%% -----------------------------------------------------------------
fs=44100; 
N=2^14;disp(['Total duration [s] ' num2str(N/fs)])
time=([1:N]'-1)/fs; % make time positive column vector starting at zero
signal=zeros(N,1);
signal(1)=1;%unit impulse

%---------Analog A-weighting filter according to IEC/CD 1672---------
% output=Aweighting_filter(signal,fs);

%---------Analog A-weighting filter according to IEC/CD 1672---------
f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
%f4 = 12194.217;%for infinite fs
f4 = 14100;%for fs=96kHz, makes magnitude -9.5dB at 20kHz
A1000 = 1.9997;
NUM = [(2*pi*f4)^2*10^(A1000/20) 0 0 0 0 ];% s^8 character
DEN = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
DEN = conv(conv(DEN,[1 2*pi*f3]),[1 2*pi*f2]);
% Bilinear transformation of analog design to get the digital filter. 
[b,a] = bilinear(NUM,DEN,fs);
%-----------------------------------------
output2=filter(b,a,signal);

figure(300)
plot(time,output,'b');
grid on;
axis([0 .003 -.1 .1]);
title('A-wtd impulse signal')
xlabel('Time [s]')
%---------------frequency domain----------------------
f=([1:N/2+1]'-1)*fs/N;
OUTPUT=fft(output);
OUTPUT2=fft(output2);

figure(400);
semilogx(f,20*log10(abs(OUTPUT2(1:floor(N/2+1)))),'b');
grid on;hold on;
semilogx(f,20*log10(abs(OUTPUT(1:floor(N/2+1)))),'g');
semilogx(fr,dBAtable,'r');
pk=10*ceil(max(20*log10(abs(OUTPUT)))/10);
axis([fs/N,fs/2,pk-60,pk])
legend('A-wtg 44.1','A-wtg 44.1 f5 f6','A-wtg table','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('bilinear magnitude response');


%% ------------------44.1 jule walker implementation-----------------
fs=44100;
Wn=2*fr/fs;
Amag=10.^(dBAtable/20);
[b,a]=yulewalk(10,Wn,Amag);% not in MATLAB documentation?
output4=filter(b,a,signal);
[d,c]=butter(1,2*240/fs,'high');% this corrects DC-LF with highpass
output=filter(d,c,output4);
[k,h]=butter(1,2*150/fs,'high');% this corrects DC-LF with highpass
output=filter(k,h,output);
[m,l]=butter(1,2*40/fs,'high');% this corrects DC-LF with highpass
%output4=filter(m,l,output);
%------------------compactify into single filter----------------
bk=conv(b,d);ah=conv(a,c);% YW+HP1
bm=conv(bk,k);al=conv(ah,h);% YW1+HP2
bo=conv(bm,m);an=conv(al,l);% YW2+HP3 final coefficients
outputconv=filter(bo,an,signal);
%-----------------plot time---------------
time=([1:N]'-1)/fs;
figure(500)
plot(time,outputconv,'b');
grid on;
axis([0 .003 -1 1]);
title('A-wtg julewalk 44.1k impulse signal')
xlabel('Time [s]')
%----------------------plot A-wtg table & yulewalk response------------
f=[0:N/2]*fs/N;
OUTPUT4=fft(output4);
OUTPUTconv=fft(outputconv);
%--------------------------------------------------------------------
figure(550);
semilogx(f,20*log10(abs(OUTPUT4(1:floor(N/2+1)))),'b');
grid on;hold on;
semilogx(f,20*log10(abs(OUTPUTconv(1:floor(N/2+1)))),'k');
semilogx(fr,dBAtable,'r');
pk=10*ceil(max(20*log10(abs(OUTPUT)))/10);
axis([fs/N,fs/2,ylim])
legend('A yulewalk','A yulewalkconv','A-wtg table','Location','Best');%,'ITUMAG','ITUFILT','ITUFILT2','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('44.1 dBA magnitude response');
%% ------------------96kHz jule walker implementation-----------------
%% -------------------------A-wtg table-----------------------------------
fr=[0 6.3 10 16 25 40 63 100 160 250 400 630 1000 1600 2500 4000 6300 10000 16000 20000 25000 35000 48000];
dBAtable=[-inf -85.4 -70.4 -56.7 -44.7 -34.6 -26.2 -19.1 -13.4 -8.6 -4.8 -1.9 0 1.0 1.3 1.0 -0.1 -2.5 -6.6 -9.3 -15 -25 -40];
fs=96000;
Wn=2*fr/fs;
Amag=10.^(dBAtable/20);
[b,a]=yulewalk(12,Wn,Amag);
output4=filter(b,a,signal);
[d,c]=butter(1,2*280/fs,'high');% this corrects DC-LF with highpass
output=filter(d,c,output4);
[k,h]=butter(1,2*200/fs,'high');% this corrects DC-LF with highpass
output=filter(k,h,output);
[m,l]=butter(1,2*40/fs,'high');% this corrects DC-LF with highpass
%output4=filter(m,l,output);
%------------------compactify into single filter----------------
bk=conv(b,d);ah=conv(a,c);% YW+HP1
bm=conv(bk,k);al=conv(ah,h);% YW1+HP2
bo=conv(bm,m);an=conv(al,l);% YW2+HP3 final coefficients
outputconv=filter(bo,an,signal);
%-----------------plot time---------------
time=([1:N]'-1)/fs;
figure(600)
plot(time,outputconv,'b');
grid on;
axis([0 .003 -1 1]);
title('A-wtg julewalk 96k impulse signal')
xlabel('Time [s]')
%----------------------plot A-wtg table & yulewalk response------------
f=[0:N/2]*fs/N;
OUTPUT4=fft(output4);
OUTPUTconv=fft(outputconv);
%--------------------------------------------------------------------
figure(650);
semilogx(f,20*log10(abs(OUTPUT4(1:floor(N/2+1)))),'b');
grid on;hold on;
semilogx(f,20*log10(abs(OUTPUTconv(1:floor(N/2+1)))),'k');
semilogx(fr,dBAtable,'r');
pk=10*ceil(max(20*log10(abs(OUTPUT)))/10);
axis([fs/N,fs/2,ylim])
legend('A yulewalk','A yulewalkconv','A-wtg table','Location','Best');%,'ITUMAG','ITUFILT','ITUFILT2','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('96k dBA magnitude response');
%disp('-----------------------------finished------------------------')


