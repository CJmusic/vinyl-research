clear all; clc;close all;
disp('----------------start of program--------------------')
set(0,'DefaultLineLinewidth',1.5)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
%


[refLockout, fs]= audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\lockoutgroove\lockout.wav');
[dataLockout, fs] = audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\lockoutgroove\lockout2.wav');

refLockout = refLockout(:,1);
dataLockout = dataLockout(:,1);

[acor_L,lags_L] = xcorr(refLockout,dataLockout);
[M_L,I_L] = max(abs(acor_L));
lagdiff_L = lags_L(I_L);
lagdiff = lagdiff_L;

disp(strcat('lagdiff ...', num2str(lagdiff)))

timeref = (0:length(ref)-1)/fs;
timedata = (0:length(data)-1)/fs  + lagdiff/fs;

figure(1)
plot(lags_L,acor_L)
figure(2)
plot(timeref, refLockout)
plot(timedata, dataLockout)