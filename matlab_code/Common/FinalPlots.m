clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/Code/vinyl-research/matlab_code/')


A0000B0000 = readtable('/Users/cz/Code/vinyl-research/matlab_code/A0000B0000/Tbl4.csv');
A0137B0137 = readtable('/Users/cz/Code/vinyl-research/matlab_code/A0137B0137/Tbl4.csv');

head(A0000B0000)
head(A0137B0137)




function plot_histogram(Tbl,plotnum,track,measurement)

end


function plot_scatter2(plotnum, x1, y1, x2, y2, titlestring)
    % figure(plotnum);  
    fig = figure('visible','off');
    scatter(x1,y1,'ko');
    grid on; hold on;
    scatter(x2,y2,'kx');
    legend('left channel', 'right channel');
    title(titlestring);
    plotfile = strcat('plots/',titlestring);
    saveas(fig, plotfile);
end

function plot_scatter(plotnum, x1, y1, titlestring)
    % figure(plotnum);  
    fig = figure('visible','off');
    scatter(x1,y1,'ko');
    grid on;
    title(titlestring);
    plotfile = strcat('plots/',titlestring);
    saveas(fig, plotfile);
end

function plot_scatteravg(plotnum, x1, y1, titlestring)
    % figure(plotnum);  
    fig = figure('visible','off');
    scatter(x1,y1,'ko');
    grid on; hold on;;
    title(titlestring);
    plotfile = strcat('plots/',titlestring);
    saveas(fig, plotfile);
end