clear all;close all;clc
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code')
addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')

% folder = ('/Volumes/AUDIOBANK/audio_files/A0000B0000/')
% plot_folder = '/Volumes/AUDIOBANK/audio_files/A0000B0000/plots/';

folder = ('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
plot_folder = '/Volumes/AUDIOBANK/audio_files/A0137B0137/plots/';


files = dir(fullfile(folder,'*.wav'))
plotfiles = dir(fullfile(plot_folder))
cellplotfiles = struct2cell(plotfiles)
plotnames = cellplotfiles{1,:}

for i = (1:length(files)) %%loop through records
    filename = files(i).name;
    files(i);
    recordnum = filename(1:end-4)

    plotfiles
    strcat(recordnum,'.png')

    if any(strcmp(plotnames, strcat(recordnum,'.png')))
        disp('...already processed')
        continue
    end
    % filename = filename(3:end);
    disp(['opening file...:', filename])
    [data,fs] = audioread(filename);
    figure(1)
    subplot(2,1,1)
    time = (1:length(data))/fs;
    plot(time,data(:,1),'k')
    title('left channel')
    subplot(2,1,2)
    plot(time,data(:,2),'k')
    title('right channel')
    plotname = strcat(plot_folder, recordnum,'.png')
    saveas(figure(1),plotname)

end