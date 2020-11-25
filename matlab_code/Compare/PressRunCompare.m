clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

c = gray(20);


% Load in A0000B0000 table/audio table
% load in A0137B0137 table/audio table

% do histograms for RMS of each

% A0000B0000 = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000-AudioTableOct26.csv');
A0000B0000 = readtable('/Users/cz/Code/vinyl-research/matlab_code/A0000B0000/Tbl4.csv');
A0137B0137 = readtable('/Users/cz/Code/vinyl-research/matlab_code/A0137B0137/Tbl4.csv');
% A0137B0137 = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137-AudioTableOct21.csv');


head(A0000B0000)
head(A0137B0137)


plot_stereohistogram(A0000B0000, 'quiet','a','RMS_L','RMS_R','First pressing RMS levels in the quiet tracks side a','A0000B0000quietRMSa')
plot_stereohistogram(A0137B0137, 'quiet','a','RMS_L','RMS_R','Second pressing RMS levels in the quiet tracks side a','A0137B0137quietRMSa')

plot_stereohistogram(A0000B0000, 'quiet','b','RMS_L','RMS_R','First pressing RMS levels in the quiet tracks side a','A0000B0000quietRMSa')
plot_stereohistogram(A0137B0137, 'quiet','b','RMS_L','RMS_R','Second pressing RMS levels in the quiet tracks side a','A0137B0137quietRMSa')



function plot_stereohistogram(Tbl, trackname, side, x1, x2, titlestring, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx1 = find(ismember(cols, x1));
    colx2 = find(ismember(cols, x2));
    data_L = table2array(Tbl(:,colx1));
    data_R = table2array(Tbl(:,colx2));


    statsL = datastats(data_L);
    statsR = datastats(data_R);    
    
    lower_binL = statsL.mean - 10;
    lower_binR = statsR.mean - 10;
    upper_binL = statsL.mean + 10;
    upper_binR = statsR.mean + 10;


    fig = figure('Visible', 'off');
    histogram(data_L, 50,'BinLimits',[lower_binL,upper_binL], 'facecolor',[0.6 0.6 0.6])
    hold on; grid on;
    histogram(data_R,50,'BinLimits',[lower_binL,upper_binL], 'facecolor',[0.3 0.3 0.3])
    title(titlestring)
    legend('left channel', 'right channel')
    dim = [0.2 0.5 0.3 0.3];
    str = {strcat('left mean :',num2str(statsL.mean)),strcat('left std :',num2str(statsL.std)),strcat('right mean :',num2str(statsR.mean)),strcat('right std :',num2str(statsR.std))};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor', 'white');
    plotname = strcat('plots/',filename,'.png');
    saveas(fig, plotname);
end