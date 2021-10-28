% Record Surface Noise analysis
% using correlation & coherence between different recordings same groove
% John Vanderkooy
% Feb. 2019
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

addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')

tracks1 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0000B0000/040318_A0000B0000r005a1555.796.wav')
tracks2 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/010a1600.264.wav')

trackname = 'transition'

data1 = tracks1(trackname);
data2 = tracks2(trackname);

% clicks_timestamps = [11.9574, 11.7464, 11.8166,11.5757, 11.9574]; %transition
% clicks_timestamps = [10.0672, 10.0017, 9.99723, 9.91687, 10.1332];  %quiet
fs = 96000; 



time = linspace(0,(length(data1)-1)/fs,length(data1));
rotation_speed = 33.33333;%45;
n_sam = fs*1.8;
num_segs = (floor(length(data1)/n_sam));
time_seg = time(1:n_sam);
seg_array = []; 

for ng = 1:num_segs
    seg_array1(:,:,ng) = data1(1+(ng-1)*n_sam:ng*n_sam,:);
    seg_array2(:,:,ng) = data2(1+(ng-1)*n_sam:ng*n_sam,:);
end

for ng = 4:5%num_segs-1
    %%% coherences in the left channel 
    disp('NEW LOOP~~~~~~~~~~~~~~~~~~~~~~~~~')
    %%********* PRINT THE LAGDIFF BETWEEN GROOVES ************
    if ng > 1;
        groovediff1 = audio_corrlineup(seg_array1(:,1,ng), seg_array1(:,1,ng-1));
        seg_array1(:,2,ng-1) = circshift(seg_array1(:,2,ng-1), groovediff1);

        groovediff2 = audio_corrlineup(seg_array2(:,1,ng), seg_array2(:,1,ng-1));
        seg_array2(:,2,ng-1) = circshift(seg_array2(:,2,ng-1), groovediff2);
    end

    %%%%%
    % coh_LR = audio_mscohere(seg_array(:,1,ng), seg_array(:,2,ng),fs); 

    [coh_nextL1, freq_coh] = audio_mscohere(seg_array1(:,1,ng), seg_array1(:,1,ng+1), fs); % previous left channel with next left channel
    [coh_firstL1, ~] = audio_mscohere(seg_array1(:,1,1), seg_array1(:,1,ng+1), fs); % first left channel with next left channel
    %%% coherences in the right channel    
    [coh_nextR1, ~] = audio_mscohere(seg_array1(:,2,ng), seg_array1(:,2,ng+1), fs);% previous right channel with next right channel
    [coh_firstR1, ~] = audio_mscohere(seg_array1(:,2,1), seg_array1(:,2,ng+1), fs);% first right channel with next right channel
    [coh_nextLR1, ~] = audio_mscohere(seg_array1(:,1,ng), seg_array1(:,2,ng+1), fs); % previous left channel with next right channel
    [coh_nextRL1, ~] = audio_mscohere(seg_array1(:,2,ng), seg_array1(:,1,ng+1), fs); % previous right channel with next left channel
    [coh_LR1, ~] = audio_mscohere(seg_array1(:,1,ng), seg_array1(:,2,ng), fs);% same groove, left and right channel

    coh_firstL1 = pwroctsmooth_singlesided(coh_firstL1, 0.33);
    coh_firstR1 = pwroctsmooth_singlesided(coh_firstR1, 0.33);
    coh_nextL1 = pwroctsmooth_singlesided(coh_nextL1, 0.33);
    coh_nextR1 = pwroctsmooth_singlesided(coh_nextR1, 0.33);
    coh_nextLR1 = pwroctsmooth_singlesided(coh_nextLR1, 0.33);
    coh_nextRL1 = pwroctsmooth_singlesided(coh_nextRL1, 0.33);
    coh_LR1 = pwroctsmooth_singlesided(coh_LR1, 0.33);


    [coh_nextL2, freq_coh] = audio_mscohere(seg_array2(:,1,ng), seg_array2(:,1,ng+1), fs); % previous left channel with next left channel
    [coh_firstL2, ~] = audio_mscohere(seg_array2(:,1,1), seg_array2(:,1,ng+1), fs); % first left channel with next left channel
    %%% coherences in the right channel    
    [coh_nextR2, ~] = audio_mscohere(seg_array2(:,2,ng), seg_array2(:,2,ng+1), fs);% previous right channel with next right channel
    [coh_firstR2, ~] = audio_mscohere(seg_array2(:,2,1), seg_array2(:,2,ng+1), fs);% first right channel with next right channel
    [coh_nextLR2, ~] = audio_mscohere(seg_array2(:,1,ng), seg_array2(:,2,ng+1), fs); % previous left channel with next right channel
    [coh_nextRL2, ~] = audio_mscohere(seg_array2(:,2,ng), seg_array2(:,1,ng+1), fs); % previous right channel with next left channel
    [coh_LR2, ~] = audio_mscohere(seg_array2(:,1,ng), seg_array2(:,2,ng), fs);% same groove, left and right channel

    coh_firstL2 = pwroctsmooth_singlesided(coh_firstL2, 0.33);
    coh_firstR2 = pwroctsmooth_singlesided(coh_firstR2, 0.33);
    coh_nextL2 = pwroctsmooth_singlesided(coh_nextL2, 0.33);
    coh_nextR2 = pwroctsmooth_singlesided(coh_nextR2, 0.33);
    coh_nextLR2 = pwroctsmooth_singlesided(coh_nextLR2, 0.33);
    coh_nextRL2 = pwroctsmooth_singlesided(coh_nextRL2, 0.33);
    coh_LR2 = pwroctsmooth_singlesided(coh_LR2, 0.33);



    figure(100)
    subplot(2,1,1)
    plot(seg_array1(:,1,ng))
    hold on ; grid on;
    subplot(2,1,2)
    plot(seg_array2(:,1,ng))
    hold on; grid on;

    figure(1)
    subplot(2,1,1)
    plot(freq_coh,coh_firstL1)
    hold on; 
    plot(freq_coh,coh_firstL2)
    title('left channel')
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    legend('CAF', 'Reinee')
    subplot(2,1,2)
    plot(freq_coh,coh_firstR1)
    hold on; 
    plot(freq_coh,coh_firstR2)
    title('right channel')
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    legend('CAF', 'Reinee')
    
    figure(2)
    subplot(2,1,1)
    plot(freq_coh,coh_nextL1)
    hold on; 
    plot(freq_coh,coh_nextL2)
    title('left channel')
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    legend('CAF', 'Reinee')
    subplot(2,1,2)
    plot(freq_coh,coh_nextR1)
    hold on; 
    plot(freq_coh,coh_nextR2)
    title('right channel')
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    legend('CAF', 'Reinee')

    figure(3); hold on;
    plot(freq_coh, coh_LR1)
    plot(freq_coh, coh_LR2)
    set(gca, 'XScale', 'log');
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')
    xlabel('Frequency (Hz)')
    title(strcat('Coherence in noise between Left and Right Channel'))
    legend('CAF', 'Reinee')


    

end
