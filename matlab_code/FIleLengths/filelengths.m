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
    record2 = SeperateTracks('/AUDIOBANK/audio_files/A0137B0137/000a.wav')

    % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/styluswear/040318_A0000B0000r001a.wav');

end 
if ispc() == true
    addpath('D:\Code\vinyl-research\matlab_code\A0137B0137')
    addpath('D:\Code\vinyl-research\matlab_code\')
    addpath('D:\Code\vinyl-research\matlab_code\Common')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions\')
    addpath('')

    record1 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/003a.wav');
    record2 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/000a.wav');    
end
data1 = record1('1kHz');
data1 = data1(:,1);
time = (1:length(data1))/96000;

folder = '/Volumes/AUDIOBANK/audio_files/A0137B0137/'
files = dir(fullfile(folder,'*.wav'))

for i = (1:length(files)) %%loop through records
    filename = files(i).name;
    fig = figure(i);
    plot(time, data1);
    % plot(data1)
    hold on; grid on;
    data2 = record2('1kHz');
    time = (1:length(data1))/96000;
    data2 = data2(:,1);
    plot(time, data2);
    name = string(filename(1:4))
    saveas(fig, name, 'png')
end
% plot(data2)


% figure(2)
% data3 = record1('leadout');
% data3 = data3(:,1);
% time = (1:length(data3))/96000;
% plot(time,data3)
% hold on; grid on;
% data4 = record2('leadout');
% data4 = data4(:,1);
% time = (1:length(data4))/96000;
% plot(time,data4)

% size(data1)
% size(data2)
% figure(3)
% plot((data1 - data2))
% sum((data1 - data2).^ 2)
% correlation = sqrt(sum((data1 - data2).^ 2))/length(data1)  % >= R2016b: auto-expand
% correlation = correlation / max(correlation)  % Normalize to [0, 1]

figure(3)
data3 = record1('sweepV2');
data3 = data3(:,1);
time = (1:length(data3))/96000;
plot(time,data3)
hold on; grid on;
data4 = record2('sweepV2');
data4 = data4(:,1);
time = (1:length(data4))/96000;
plot(time,data4)
