clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


files = dir('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/antiskatetest/*.wav')

RMS_L = []
RMS_R = []
AS = []

for i = (1:length(files))
    filename = files(i).name
    file = strcat(files(i).folder,'/',filename)

    as = str2num(filename(1:3))
    AS = [AS, as]

    [data, fs] = audioread(file); 
    % data = data(5*fs:15*fs);

    rms_l = rms(data(:,1))
    rms_r = rms(data(:,2))

    RMS_L = [RMS_L, rms_l];
    RMS_R = [RMS_R, rms_r];
    plot(data)

end
% H = bar(AS,[RMS_L, RMS_R], 'LineWidth', 2)
H = bar([RMS_L, RMS_R], 'LineWidth', 2)
RMS_L
RMS_R
AS
H(1).FaceColor = [0.6 0.6 0.6];
H(2).FaceColor = [.9 .9 .9];
