clear all;close all;clc
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

if ismac() == true
    addpath('/Users/cz/Code/vinyl-research/matlab_code')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/from_John')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/A0137B0137')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
    % addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/')

    % reference = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/samerecordtest/3.wav');
    % record = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/samerecordtest/2.wav');

    % reference = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav');
    % record = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/styluswear/042820_A0000B0000r1a.wav');
    % reference = SeperateTracks('/AUDIOBANK/audio_files/A0137B0137/003a.wav')
    % disp('Calling seperatetracks')
    % reference = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav')

    % reference = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r029a.wav')


    % record = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3b.wav')

    % record = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/styluswear/040318_A0000B0000r001a.wav');

    % reference = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r029a.wav')
    % record = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav')

    reference = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/003a.wav')
    record = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/000b.wav')



end 
if ispc() == true
    addpath('D:\Code\vinyl-research\matlab_code\A0137B0137')
    addpath('D:\Code\vinyl-research\matlab_code\')
    addpath('D:\Code\vinyl-research\matlab_code\Common')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions\')
    addpath('')

    % reference = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/003a.wav');
    reference = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/-01b.wav'); offset1 = 957.35208;
    % record = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/000a.wav');    
    % record = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/012b.wav'); offset2 = 957.5021;
    record = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/-01b.wav'); 
end



refsweep = reference('sweep');
refsweep2 = reference('sweep2');
sweep = record('sweep');
sweep2 = record('sweep2');

refleadout = reference('leadout');
leadout = record('leadout');

fs = 96000;


figure(4)
[acor_L,lags_L] = xcorr(refleadout(:,1),leadout(:,1));
[M_L,I_L] = max(abs(acor_L));
lagdiff_L = lags_L(I_L);
lagdiff = lagdiff_L;
plot(lags_L/fs, acor_L, 'k')
% xlim([-1, 1])
grid on
title('cross correlation of leadouts')
saveas(figure(4),'xcorrleadouts.png')


figure(1)
[acor_L,lags_L] = xcorr(refsweep(:,1),sweep(:,1));
[M_L,I_L] = max(abs(acor_L));
lagdiff_L = lags_L(I_L);
lagdiff = lagdiff_L;
plot(lags_L/fs, acor_L, 'k')
% xlim([-1, 1])
grid on
title('cross correlation of sweeps')
% plotname = strcat(plot_folder, recordnum,'sweep.png')
saveas(figure(1),'xcorrsweeps.png')

figure(2)
[acor_L,lags_L] = xcorr(refsweep2(:,1),sweep2(:,1));
[M_L,I_L] = max(abs(acor_L));
lagdiff_L = lags_L(I_L);
lagdiff = lagdiff_L;
plot(lags_L/fs, acor_L, 'k')
% xlim([-1, 1])
grid on
title('cross correlation of sweeps 2')
% plotname = strcat(plot_folder, recordnum,'sweep2.png')
saveas(figure(2),'xcorrsweeps2.png')

figure(3)
time = (1:length(refleadout))/96000;
plot(time,refleadout(:,1))
hold on; grid on;
time = (1:length(leadout))/96000;
plot(time,leadout(:,1))
title('leadouts')
% plotname = strcat(plot_folder, recordnum,'leadouts.png')
saveas(figure(3), 'leadouts.png')

figure(5)
time = (1:length(refsweep))/96000;
plot(time,refsweep(:,1))
hold on; grid on;
time = (1:length(sweep))/96000;
plot(time,sweep(:,1))
title('sweeps')
% plotname = strcat(plot_folder, recordnum,'leadouts.png')
saveas(figure(5), 'sweep.png')


% data1 = reference('1kHz');
% data1 = data1(1:length(data1)-1,1);
% time = (1:length(data1))/96000;

% data1 = data1/max(data1);
% plot(time, data1);

% % plot(data1)
% hold on; grid on;
% data2 = record('1kHz');
% time = (1:length(data1))/96000;
% data2 = data2(1:length(data1),1);
% data2 = data2/max(data2);
% size(data1)
% size(data2)
% plot(time, data2);
% % plot(data2)


% figure(2)
% data3 = reference('leadout');
% data3 = data3(:,1);
% time = (1:length(data3))/96000;
% plot(time,data3)
% hold on; grid on;
% data4 = record('leadout');
% data4 = data4(:,1);
% time = (1:length(data4))/96000;
% plot(time,data4)
% title('leadouts')
% saveas(figure(2), 'leadouts.png')

% figure(3)
% data3 = reference('sweepV2');
% data3 = data3(:,1);
% time = (1:length(data3))/96000;
% plot(time,data3)
% hold on; grid on;
% data4 = record('sweepV2');
% data4 = data4(:,1);
% time = (1:length(data4))/96000;
% plot(time,data4)
% title('sweeps')
% saveas(figure(3), 'sweeps.png')

% fs = 96000;
% datastart = 58*fs;
% dataend = 60*fs;
% data1 = data1(datastart:dataend);

% figure(4)
% [b,a]=butter(2,[700*2/fs 1400*2/fs]);
% data1=filtfilt(b,a,data1);
% signal = hilbert(data1);
% time = (1:length(data1))/96000;
% plot(time,abs(signal))
% N = length(data1)
% startsum=round(0.01*N);stopsum=round(0.25*N);
% amplitde=sum(abs(signal(startsum:stopsum)))/(stopsum-startsum+1)
% grid on;

% figure(5)
% [acor_L,lags_L] = xcorr(data3,data4);
% [M_L,I_L] = max(abs(acor_L));
% lagdiff_L = lags_L(I_L);
% lagdiff = lagdiff_L;
% plot(lags_L/fs, acor_L, 'k')
% xlim([-5, 5])
% grid on
% title('cross correlation of sweeps')
% saveas(figure(5),'xcorr.png')