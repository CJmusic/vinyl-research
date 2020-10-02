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
    addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')

    % record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/samerecordtest/3.wav');
    % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/samerecordtest/2.wav');

    % record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav');
    % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/styluswear/042820_A0000B0000r1a.wav');
    reference = SeperateTracks('/AUDIOBANK/audio_files/A0137B0137/003a.wav')
    % record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/ableton/r27a-ableton1024a.wav')
    folder = '/Volumes/AUDIOBANK/audio_files/A0137B0137/'
    % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/styluswear/040318_A0000B0000r001a.wav');

end 
if ispc() == true
    addpath('D:\Code\vinyl-research\matlab_code\A0137B0137')
    addpath('D:\Code\vinyl-research\matlab_code\')
    addpath('D:\Code\vinyl-research\matlab_code\Common')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions\')
    addpath('')

    reference = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/003a.wav');
    % record2 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/000a.wav');    
end

% data1 = record1('1kHz');
% data1 = data1(1:length(data1)-1,1);
% time = (1:length(data1))/96000;

folder = '/Volumes/AUDIOBANK/audio_files/A0137B0137/'
files = dir(fullfile(folder,'*.wav'))

plot_folder = '/Users/cz/Code/vinyl-research/matlab_code/FileLengths/plots/';
fs = 96000;
for i = (1:length(files)) %%loop through records
    filename = files(i).name
    record = SeperateTracks(filename);
    recordnum = filename(1:end-4)

    refsweep = reference('sweep');
    refsweep2 = reference('sweep2');
    sweep = record('sweep');
    sweep2 = record('sweep2');

    refleadout = reference('leadout');
    leadout = record('leadout');


    figure(1)
    [acor_L,lags_L] = xcorr(refsweep(:,1),sweep(:,1));
    [M_L,I_L] = max(abs(acor_L));
    lagdiff_L = lags_L(I_L);
    lagdiff = lagdiff_L;
    plot(lags_L/fs, acor_L, 'k')
    xlim([-1, 1])
    grid on
    title('cross correlation of sweeps')
    plotname = strcat(plot_folder, recordnum,'sweep.png')
    saveas(figure(1),plotname)
    
    figure(2)
    [acor_L,lags_L] = xcorr(refsweep2(:,1),sweep2(:,1));
    [M_L,I_L] = max(abs(acor_L));
    lagdiff_L = lags_L(I_L);
    lagdiff = lagdiff_L;
    plot(lags_L/fs, acor_L, 'k')
    xlim([-1, 1])
    grid on
    title('cross correlation of sweeps 2')
    plotname = strcat(plot_folder, recordnum,'sweep2.png')
    saveas(figure(2),plotname)

    figure(3)
    time = (1:length(refleadout))/96000;
    plot(time,refleadout(:,1))
    hold on; grid on;
    time = (1:length(leadout))/96000;
    plot(time,leadout(:,1))
    title('leadouts')
    plotname = strcat(plot_folder, recordnum,'leadouts.png')
    saveas(figure(3), plotname)

    % %Plot the two sweeps on top of each other
    % fig = figure(1)
    % data3 = record1('sweepV2');
    % data3 = data3(:,1);
    % time = (1:length(data3))/96000;
    % % plot(time,data3)
    % % hold on; grid on;
    % data4 = record2('sweepV2');
    % data4 = data4(:,1);
    % time = (1:length(data4))/96000;
    % plot(time,data4)

    % %-----------------crosscorrelations----------
    % nlags=100000;
    % fs = 96000;
    % [C1,tlag]=xcorr(data3,data4,nlags);
    % % [C2,tlag]=xcorr(dut2,ref2,nlags);
    % %----------plot amplitude------------
    % figure(20)
    % subplot(2,1,1)

    % data3 = data3(:,1);
    % time = (1:length(data3))/96000;
    % plot(time,data3)
    % grid on; hold on;
    % data4 = record2('sweepV2');
    % data4 = data4(:,1);
    % time = (1:length(data4))/96000;
    % plot(time,data4)
    % subplot(2,1,2)
    % plot(tlag/fs,C1,'b')
    % grid on; hold on;
    % legend('1st sweep','2nd sweep')
    % xlabel('time [s]')
    % ylabel('crosscorrelation')
    % axis([xlim ylim])
    % title('very coarse comparison')

    % clf(figure(20))
    % clf(figure(1))
    % %----------plot amplitude------------
    % % figure(30)
    % % plot(tlag/fs,C1,'b')
    % grid on;hold on
    % legend('1st sweep','2nd sweep')
    % xlabel('time [s]')
    % ylabel('crosscorrelation')
    % axis([xlim/10 ylim])
    % title('positive offset means dut missing samples wrt reference')
    % %----------plot amplitude------------
    % figure(40)
    % plot(tlag/fs,C1,'b')
    % grid on;hold on
    % legend('1st sweep','2nd sweep')
    % xlabel('time [s]')
    % ylabel('crosscorrelation')
    % axis([xlim/100 ylim])
    % title('very fine comparison')
    % disp('-------------end of test------------------')
    % figname = strcat(folder, filename(2:end-4),'.png')
    % saveas(fig,figname)
    % clf(fig)
    clf(figure(1))
    clf(figure(2))
    clf(figure(3))
end


