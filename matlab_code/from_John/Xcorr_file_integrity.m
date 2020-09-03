% Record file integrity analysis using xcorr on sweep
% July 2020 John Vanderkooy
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
%filename='27b_declick.wav';t1=956.8037;% example of bad file
%filename='r29a-Technics-full.wav';t1=957.7999;
%filename='r30a-Technics-full.wav';t1=958.73495;
%filename='r30b-Technics-full.wav';t1=952.7588;

% filename='r31a-Technics-full.wav';t1=955.48471;
filename='D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\040318_A0000B0000r031a.wav';t1=954.1865;

%filename='r31b-Technics-full.wav';t1=958.9805;
%filename='r104b-1-declicked.wav';t1=958.4144;
%filename='minpucksize1b.wav';t1=957.352;
filename='D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0137B0137\010b.wav';t1=957.352;

%filename='minpucksizem6b.wav';t1=957.5021;
%filename='042820_late_test_file.wav';t1=957.2155;%*********reference file
%filename='042820_late_test_file_declicked.wav';t1=957.2155;
%filename='040318_early_test_file.wav';t1=950.7134;% bad file

reffilename='D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\styluswear\042820_A0000B0000r1a.wav';t1ref=957.2155;%*********reference file
% reffilename='042820_late_test_file.wav';t1ref=957.2155;%*********reference file

[sigref,fs]=audioread(reffilename);
[signal,fs]=audioread(filename);

refsweepstart=137;%duration =35sec
refsweepstart2=655;%inner sweep
%----------find indices for analysis-----------
fileend=round(t1*fs);%fiducial lockout groove for file under test
refend=round(t1ref*fs);%fiducial lockout groove for reference file
refsweepstart=round(655*fs);%137 or 655 sec
offset=refend-refsweepstart;
sigsweepstart=fileend-offset;
%----------select reference and signal portions----------
lr=1;
sigref=sigref(refsweepstart:refsweepstart+round(35*fs),lr);
signal=signal(sigsweepstart:sigsweepstart+round(35*fs),lr);

N=length(signal)%length in samples of audio file
t=(0:N-1)/fs;
%------------------
figure(10)
plot(t,sigref)
grid on
xlabel('time [s]')
ylabel('ref sweep')
title('reference')
%------------------
figure(15)
plot(t,signal)
grid on
xlabel('time [s]')
ylabel('sig sweep')
title('signal')
%-----------------crosscorrelation----------
[C,tlag]=xcorr(sigref,signal,10000);
%----------plot amplitude------------
figure(20)
plot(tlag/fs,C)
grid on
xlabel('time [s]')
ylabel('crosscorrelation')
disp('-------------end of test------------------')




