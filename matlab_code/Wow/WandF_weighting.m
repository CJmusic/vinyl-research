%W&F-weighting digital filter
% this program optimizes the parameters of an s^3 HP x 1st LP

clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
try
    pkg load signal; %for Octave
catch
end

%--------------------------------------
fs=9600;%decimated down from 96000
N=2^20;disp(['Total duration [s] ' num2str(N/fs)])
time=([1:N]'-1)/fs; % make time positive column vector starting at zero
signal=zeros(N,1);
signal(10)=1;%unit impulse


%% -------------------------WF-wtg table-----------------------------------
fr=[0.1 0.19 0.43 0.77 1.0 2.0 5.0 10.0 20.0 50.0 165 1000];
dBWFtable=[-57 -40 -20 -10 -7.25 -1.52 0 -1 -4 -10 -20 -36];
%---------------------------------------

figure(60);
semilogx(fr,dBWFtable,'k');
grid on;hold on
%axis([fs/N,fs/2,pk-60,pk])
% legend(' W&F-wtg table','Location','Best')
xlabel('Frequency [Hz]')
ylabel('dB')
% title('wow and flutter weighting');
saveas(figure(60),'wowfiltfreqresponse.png')


f1 = 15.0;%HF rolloff
f2 = 0.65;%LF rollup
f3 = 0.9;%LF rollup
f4 = 1.;%LF rollup
WF4 = 0.71;%sets dB gain
X=[f1 f2 f3 f4 WF4]

%---------Analog W&F-weighting filter from filter convolution---------

NUM = X(5)*[(2*pi)^3*X(2)*X(3)*X(4) 0 0 0];% s^3 character
DEN = conv(conv(conv([1 2*pi*X(2)],[1 2*pi*X(3)]) ,[1 2*pi*X(4)]), [1 2*pi*X(1)]); 
% Bilinear transformation of analog design to get the digital filter. 
[b,a] = bilinear(NUM,DEN,fs);
%-----------------------------------------
output=filter(b,a,signal);
%-----------------get dB frequency response error----------
h=freqz(b,a,fr,fs);
hdB=20*log10(abs(h));
er=hdB-dBWFtable;


%---------------frequency domain----------------------
f=[0:N/2]'*fs/N;
OUTPUT=fft(output);

figure(400);
semilogx(f,20*log10(abs(OUTPUT(1:floor(N/2+1)))),'b');
grid on;hold on;
semilogx(fr,dBWFtable,'r');
pk=ceil(max(dBWFtable));
axis([10*fs/N,fs/9.6,pk-60,pk+20])
legend('WF-bilinear','WF-wtg table','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('weighting magnitude responses');

disp('-----------------------------finished------------------------')


