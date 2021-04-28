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

addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')


[data, fs] = audioread('/Users/cz/Documents/for_John2/028b.wav');
[data2, fs] = audioread('/Users/cz/Documents/for_John2/032b.wav');


t_dropout = 847.477;
t_dropout = 847.477 - 45/fs;

data_dropout = data((t_dropout - 0.001)*fs:(t_dropout + 0.001)*fs,:);

fig = figure(1);
% time = (1:length(data_dropout))/fs;
stem(data_dropout(:,2),'k');
ylabel('Signal Level')
xlabel('Samples')
title('Plot of a sample dropout')
% hold on; grid on;
saveas(fig, 'plots/sampledropout.png')


data = data(:,:)/rms(data(28*fs:30*fs));
data2 = data2(:,:)/rms(data2(28*fs:30*fs));

size(data)
size(data2)

fig = figure(2);
ts = 831.420; tf = 853.146;
time = (1:(tf-ts)*fs+1)/fs;
subplot(2,1,1)
plot(time,data(ts*fs:tf*fs),'k')
ylabel('Signal Level')
xlabel('Time [s]')
ylim([-0.5,0.5])
xlim([0,tf-ts])
title('Record 028 side b 1kHzR track right channel')
subplot(2,1,2)
plot(time,data2(ts*fs:tf*fs),'k')
ylabel('Signal Level')
xlabel('Time [s]')
ylim([-0.5,0.5])
xlim([0,tf-ts])
title('Record 032 side b 1kHzR track right channel')
saveas(fig, 'plots/tracklength.png')
