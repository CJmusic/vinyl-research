% Record RMS Measurement
% using rms, filter, declick
% John Vanderkooy
% Oct 2020
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

filename='042820_late_test_file_1kHz.wav';ts=35;tf=45;
[rev,fs]=audioread(filename);

w1=2*707/fs;w2=2*1404/fs;%bandpass for 1kHz
[b,a]=butter(4,[w1 w2]);
revf=filter(b,a,rev);

filename='042820_late_test_file_declick_1kHz.wav';%ts=20;tf=30;
[revdeclick,fs]=audioread(filename);

lr=2; %1=left, 2=right
disp(['left_right: ' num2str(lr)])
Nt=length(rev);
disp(['fs: ' num2str(fs) '  N_total: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector

%-------------------plot all data------------
figure(20)
plot(t,rev(1:Nt,lr),'b') 
grid on;hold on
plot(t,revf(1:Nt,lr),'r')
plot(t,revdeclick(1:Nt,lr),'g') 
xlabel('time[s]')
legend('raw','filt','declicked')
title('untrimmed file data')
%----------select useful portion------------
ns=round(ts*fs)+1;nf=round(tf*fs);
rev=rev(ns:nf,lr);% after this the lr dimension is gone
revf=revf(ns:nf,lr);% after this the lr dimension is gone
revdeclick=revdeclick(ns:nf,lr);% after this the lr dimension is gone
Nt=length(rev);
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
disp(['N_analyze: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
%------------------------------
figure(30)
plot(t,rev(1:Nt),'b') 
grid on;hold on
plot(t,revf(1:Nt),'r') 
plot(t,revdeclick(1:Nt),'g')
xlabel('time[s]')
legend('raw','filt','declicked')
title('trimmed file data')
%-----------------rms determinations------------
rms=rms_response(rev)
rmsf=rms_response(revf)
rmsdeclick=rms_response(revdeclick)
%-----------------freq plots------------
f=[0:Nt/2]*fs/Nt;

REV=abs(fft(rev));
REVF=abs(fft(revf));
REVdeclick=abs(fft(revdeclick));

figure(40)
plot(f,20*log10(REV(1:floor(Nt/2+1))),'b');
grid on;hold on
plot(f,20*log10(REVF(1:floor(Nt/2+1))),'r');
plot(f,20*log10(REVdeclick(1:floor(Nt/2+1))),'g');
xlabel('freq[Hz]')
ylabel('Power Spectrum [dB]')
axis([0 2500 ylim])
legend('raw','filt','declick');
title('PSD')



disp('-------------------finished--------------------') 
