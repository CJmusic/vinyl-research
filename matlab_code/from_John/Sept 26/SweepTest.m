% Test record sweep track
% the lockout groove has clicks to sync the files at the end
% Sept. 2019 John Vanderkooy
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

%large files take time to load

%----------------shift data into sync at end-----------
fs=96000;
% test positions are read from Audacity, reference track

% filename='r104b-1.wav';t_pip=958.4143;
% filename='r104b-2.wav';t_pip=959.3487;
% filename='r104b-3.wav';t_pip=958.6623;
% filename='r104b-4.wav';t_pip=957.1195;
%filename='r104b-5.wav';t_pip=957.7983;ts=137;tf=172; %first sweep
% filename='r104b-5.wav';t_pip=957.7983;ts=655;tf=691; %third sweep
%filename='r104b-5.wav';t_pip=957.7983;ts=397;tf=432; %first vertical
filename='r104b-5.wav';t_pip=957.7983;ts=914;tf=949.5; %final vertical

[sig,fs]=audioread(filename);
Nt=length(sig)
n_pip=round(t_pip*fs);
n_start=round(ts*fs);
n_finish=round(tf*fs);
Ns=n_finish-n_start+1
sig1=sig(n_start:n_start+Ns-1,:);

clear sig;% hardly necessary now with 6GB of ram

t=[0:Ns-1]'/fs;%column vector

%-------------------plot all data------------
figure(10)
plot(t,sig1(:,1),'b') 
grid on;
axis([xlim ylim])
xlabel('time[s]')
legend('sweep file')
title('raw data')

%-----------------------data plots------------
% figure(20)
% plot(t,sig1(:,1),'b') 
% grid on;
% axis([0 0.0025*Nt/fs ylim])
% xlabel('time[s]')
% legend('file1','file2')
% title('zoom of files')

%-----------------spectra------------------
SIG1=fft(sig1(:,1));
SIG1=boxsmooth(abs(SIG1),100);
f=[0:Ns/2]'*fs/Ns;
SIG1flat=SIG1(1:floor(Ns/2+1)).*sqrt(f/fs);
%SIG1flat=SIG1;
figure(30)
semilogx(f,20*log10(abs(SIG1flat(1:floor(Ns/2+1)))),'b');
grid on;
axis([1 fs/2 ylim])
xlabel('freq[Hz]');
ylabel('[dB]');
title('Spectrum sig1')
disp('-------------------finished--------------------')
%filename='lockout-acer.wav';t_pip=599.146;% this has bad timings (dropouts)
%filename='r9technics_declicked.wav';t_pip=960.8963;% comp file
%filename='071819-A0000B000r24a.wav';t_pip=956.1219;% comp file
%filename='27b_declick.wav';t_pip=956.803;% slightly off reference file%%%%%%%%%%%%%%
%filename='28b_declick.wav';t_pip=954.88387;% file to compare with reference
%filename='r28a-mac.wav';t_pip=959.3804;% fi
%filename='r28a-acer.wav';t_pip=952.9971;% fi

