% Record file integrity analysis
% June 2020 John Vanderkooy
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
%filename='27b_declick.wav';t1=956.8037;%example of bad file
%filename='r29a-Technics-full.wav';t1=957.7999;
%filename='r30a-Technics-full.wav';t1=958.73495;
%filename='r30b-Technics-full.wav';t1=952.7588;
%filename='r31a-Technics-full.wav';t1=955.48471;
%filename='r31b-Technics-full.wav';t1=958.9805;
%filename='r104b-1-declicked.wav';t1=958.4144;
%filename='minpucksize1b.wav';t1=957.352;
%filename='minpucksize1b.wav';t1=957.352;
%filename='minpucksize1b.wav';t1=957.352;
%filename='minpucksizem6b.wav';t1=957.5021;
%filename='042820_late_test_file.wav';t1=957.2155;%*********reference file
filename='042820_late_test_file_declicked.wav';t1=957.2155;
%filename='040318_early_test_file.wav';t1=950.7134;%subtract 3.5 sec extra

[signalt,fs]=audioread(filename);

%----------find indices for analysis-----------
fileend=round(t1*fs);%fiducial lockout groove for file under test
refend=round(957.2155*fs);%fiducial lockout groove for reference file
offset=fileend-refend;

nflat1kHz=offset+round(72*fs);
nstopslope=offset+round(74.0*fs);
%----------select signal portions----------
lr=1;
signal=signalt(nflat1kHz+1:nstopslope,lr);
N=length(signal)%length in samples of audio file
t=(0:N-1)/fs;
clear signalt;
%------------------
figure(10)
plot(t,signal)
grid on
xlabel('time [s]')
ylabel('zoomed signal')
%------------bandpass filter signal-------------
[b,a]=butter(2,[700*2/fs 1400*2/fs]);
signal=filtfilt(b,a,signal);
%------------------
figure(15)
plot(t,signal)
grid on
xlabel('time [s]')
ylabel('filtered signal')
%-----------------analytic signal----------
signal=hilbert(signal);
startsum=round(0.01*N);stopsum=round(0.25*N);
amplitde=sum(abs(signal(startsum:stopsum)))/(stopsum-startsum+1)
%----------plot amplitude------------
figure(20)
plot(t,abs(signal))
grid on
xlabel('time [s]')
ylabel('signal')
disp('-------------end of test------------------')




