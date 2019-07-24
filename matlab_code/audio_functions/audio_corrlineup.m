
%~ Loop through all the files and line them up with the reference

%{ LINEUP %}
% This function needs to take two arrays, and return two arrays that are lined up with one another. 
% It should: 
%       -  
%
%This function doesn't like stereo files as they come so I'll need to clean them up, then realign and redo it
%
%
%
%

clc; close all;
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/');

% audio_bin1 = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r30a-sin.wav';
% audio_bin2 = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r27a-sin.wav';

audio_bin1 = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/leadouts/03141_A0000B0000r30a.wav';
audio_bin2 = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/leadouts/03141_A0000B0000r28a.wav';

[data1, fs] = audioread(audio_bin1);
[data2, fs] = audioread(audio_bin2);

% lagdiff = audio_corrlineuptest(data1, data2, fs)
disp('LEFT THEN RIGHT')
lagdiff = audio_corrlineuptest(data1(:,1), data2(:,1), fs)
lagdiff = audio_corrlineuptest(data1(:,2), data2(:,2), fs)

data1corr = circshift(data1, lagdiff);

time1 = (0:(length(data1)-1))/fs;
time2 = (0:(length(data2)-1))/fs;


size(data1)
size(data2)
size(time1)
size(time2)


figure(1)
hold on; grid on; 
plot(time2, data2(:,1))
% plot(time2, data2(:,1), '--')
plot(time1, data1corr(:,1))


disp('datax')
size(data1(1:length(data1corr)))
size(data1corr(:,1))



figure(2)
hold on; grid on;
datax = data1(1:length(data1corr)) - data1corr(:,1).';
plot(datax)

function lagdiff = audio_corrlineuptest(data_file, data_ref, fs_file)
% function lagdiff = audio_corrlineup(data_file, data_ref)
    disp('audio_lineup function')

    % [acor_L,lags_L] = xcorr(data_file(:,1),data_ref(:,1));
    [acor_L,lags_L] = xcorr(data_ref(:,1),data_file(:,1));

    % figure(30); 
    % plot(acor_L)

    [M_L,I_L] = max(abs(acor_L));
    lagdiff_L = lags_L(I_L);
    
    % disp('lagdiff_L')
    % lagdiff_L
    
    lagdiff = lagdiff_L;

    disp ('ending audio_lineup function')
end

