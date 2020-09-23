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
    record1 = SeperateTracks('/AUDIOBANK/audio_files/A0137B0137/003a.wav')
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

    record1 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/003a.wav');
    % record2 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/000a.wav');    
end

% data1 = record1('1kHz');
% data1 = data1(1:length(data1)-1,1);
% time = (1:length(data1))/96000;

folder = '/Volumes/AUDIOBANK/audio_files/A0137B0137/'
files = dir(fullfile(folder,'*.wav'))

for i = (1:length(files)) %%loop through records
    filename = files(i).name;
    record2 = SeperateTracks(filename);
    % recordnum = filename(end-7,end-3);

    %Plot the two sweeps on top of each other
    fig = figure(1)
    data3 = record1('sweepV2');
    data3 = data3(:,1);
    time = (1:length(data3))/96000;
    % plot(time,data3)
    % hold on; grid on;
    data4 = record2('sweepV2');
    data4 = data4(:,1);
    time = (1:length(data4))/96000;
    plot(time,data4)

    %-----------------crosscorrelations----------
    nlags=100000;
    fs = 96000;
    [C1,tlag]=xcorr(data3,data4,nlags);
    % [C2,tlag]=xcorr(dut2,ref2,nlags);
    %----------plot amplitude------------
    figure(20)
    subplot(2,1,1)

    data3 = data3(:,1);
    time = (1:length(data3))/96000;
    plot(time,data3)
    grid on; hold on;
    data4 = record2('sweepV2');
    data4 = data4(:,1);
    time = (1:length(data4))/96000;
    plot(time,data4)
    subplot(2,1,2)
    plot(tlag/fs,C1,'b')
    grid on; hold on;
    legend('1st sweep','2nd sweep')
    xlabel('time [s]')
    ylabel('crosscorrelation')
    axis([xlim ylim])
    title('very coarse comparison')

    clf(figure(20))
    clf(figure(1))
    %----------plot amplitude------------
    % figure(30)
    % plot(tlag/fs,C1,'b')
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
    figname = strcat(folder, filename(2:end-4),'.png')
    saveas(fig,figname)
    clf(fig)
end


