clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/Code/vinyl-research/matlab_code/')

A0000B0000 = readtable('/Users/cz/Code/vinyl-research/matlab_code/A0000B0000/A0000B0000.csv');
A0137B0137 = readtable('/Users/cz/Code/vinyl-research/matlab_code/A0137B0137/A0137B0137.csv');

head(A0000B0000)
head(A0137B0137)

BOTH = [A0000B0000; A0137B0137];

% %~~~~~~~~~~~~ QUIET ~~~~~~~~~~%

%     plot_scatter2(A0137B0137, 'quiet', 'a', 'PressingNumber', 'clicks_L', 'clicks_R', 'A0137B0137 clicks vs pressing number side a')
%     plot_scatter2(A0137B0137, 'quiet', 'b', 'PressingNumber', 'clicks_L', 'clicks_R', 'A0137B0137 clicks vs pressing number side b')
%     plot_scatter2(A0000B0000, 'quiet', 'a', 'PressingNumber', 'clicks_L', 'clicks_R', 'A0000B0000 clicks vs pressing number side a')
%     plot_scatter2(A0000B0000, 'quiet', 'b', 'PressingNumber', 'clicks_L', 'clicks_R', 'A0000B0000 clicks vs pressing number side b')

%     plot_scatter2(A0137B0137, 'quiet', 'a', 'PressingNumber', 'commonclicksa_L', 'commonclicksa_R', 'A0137B0137 commonclicks a vs pressing number side a')
%     plot_scatter2(A0137B0137, 'quiet', 'b', 'PressingNumber', 'commonclicksa_L', 'commonclicksa_R', 'A0137B0137 commonclicks a vs pressing number side b')
%     plot_scatter2(A0000B0000, 'quiet', 'a', 'PressingNumber', 'commonclicksa_L', 'commonclicksa_R', 'A0000B0000 commonclicks a vs pressing number side a')
%     plot_scatter2(A0000B0000, 'quiet', 'b', 'PressingNumber', 'commonclicksa_L', 'commonclicksa_R', 'A0000B0000 commonclicks a vs pressing number side b')

%     plot_scatter2(A0137B0137, 'quiet', 'a', 'PressingNumber', 'RMS_L', 'RMS_R', 'A0137B0137 RMS vs pressing number side a')
%     plot_scatter2(A0137B0137, 'quiet', 'b', 'PressingNumber', 'RMS_L', 'RMS_R', 'A0137B0137 RMS vs pressing number side b')
%     plot_scatter2(A0000B0000, 'quiet', 'a', 'PressingNumber', 'RMS_L', 'RMS_R', 'A0000B0000 RMS vs pressing number side a')
%     plot_scatter2(A0000B0000, 'quiet', 'b', 'PressingNumber', 'RMS_L', 'RMS_R', 'A0000B0000 RMS vs pressing number side b')


%     plot_scatter2(A0137B0137, 'quiet', 'a', 'PressingNumber', 'A_L', 'A_R', 'A0137B0137 A weighted RMS vs pressing number side a')
%     plot_scatter2(A0137B0137, 'quiet', 'b', 'PressingNumber', 'A_L', 'A_R', 'A0137B0137 A weighted RMS vs pressing number side b')
%     plot_scatter2(A0000B0000, 'quiet', 'a', 'PressingNumber', 'A_L', 'A_R', 'A0000B0000 A weighted RMS vs pressing number side a')
%     plot_scatter2(A0000B0000, 'quiet', 'b', 'PressingNumber', 'A_L', 'A_R', 'A0000B0000 A weighted RMS vs pressing number side b')


%     plot_scatter2(A0137B0137, 'quiet', 'a', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'A0137B0137 CCIR weighted RMS vs pressing number side a')
%     plot_scatter2(A0137B0137, 'quiet', 'b', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'A0137B0137 CCIR weighted RMS vs pressing number side b')
%     plot_scatter2(A0000B0000, 'quiet', 'a', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'A0000B0000 CCIR weighted RMS vs pressing number side a')
%     plot_scatter2(A0000B0000, 'quiet', 'b', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'A0000B0000 CCIR weighted RMS vs pressing number side b')



%     a00b00RMS_La = get_stats(A0000B0000, 'quiet', 'a', 'RMS_L')
%     a00b00RMS_Rb = get_stats(A0000B0000, 'quiet', 'b', 'RMS_R')

%     a01b01RMS_La = get_stats(A0137B0137, 'quiet', 'a', 'RMS_L')
%     a01b01RMS_Rb = get_stats(A0137B0137, 'quiet', 'b', 'RMS_R')


%     a00b00clicks_La = get_stats(A0000B0000, 'quiet', 'a', 'clicks_L')
%     a00b00clicks_Rb = get_stats(A0000B0000, 'quiet', 'b', 'clicks_R')

%     a01b01clicks_La = get_stats(A0137B0137, 'quiet', 'a', 'clicks_L')
%     a01b01clicks_Rb = get_stats(A0137B0137, 'quiet', 'b', 'clicks_R')


%     a00b00commonclicks_Laa = get_stats(A0000B0000, 'quiet', 'a', 'commonclicksa_L')
%     a00b00commonclicks_Rab = get_stats(A0000B0000, 'quiet', 'b', 'commonclicksa_R')
%     a00b00commonclicks_Lba = get_stats(A0000B0000, 'quiet', 'a', 'commonclicksb_L')
%     a00b00commonclicks_Rbb = get_stats(A0000B0000, 'quiet', 'b', 'commonclicksb_R')


%     a01b01commonclicks_Laa = get_stats(A0137B0137, 'quiet', 'a', 'commonclicksa_L')
%     a01b01commonclicks_Rab = get_stats(A0137B0137, 'quiet', 'b', 'commonclicksa_R')
%     a01b01commonclicks_Lba = get_stats(A0137B0137, 'quiet', 'a', 'commonclicksb_L')
%     a01b01commonclicks_Rbb = get_stats(A0137B0137, 'quiet', 'b', 'commonclicksb_R')



% %~~~~~~~~~~~~~~~~ QUIET ENDS ~~~~~~~~~~~~~~% 

% %~~~~~~~~~~~~~ 1 kHz ~~~~~~~~~~~~~%

% plot_scatter2(A0137B0137, '1kHz', 'a', 'PressingNumber', 'RMS_L', 'RMS_R', 'A0137B0137 1kHz RMS vs pressing number side a')
% plot_scatter2(A0137B0137, '1kHz', 'b', 'PressingNumber', 'RMS_L', 'RMS_R', 'A0137B0137 1kHz RMS vs pressing number side b')
% plot_scatter2(A0000B0000, '1kHz', 'a', 'PressingNumber', 'RMS_L', 'RMS_R', 'A0000B0000 1kHz RMS vs pressing number side a')
% plot_scatter2(A0000B0000, '1kHz', 'b', 'PressingNumber', 'RMS_L', 'RMS_R', 'A0000B0000 1kHz RMS vs pressing number side b')


% plot_scatter2(A0137B0137, '1kHz', 'a', 'PressingNumber', 'A_L', 'A_R', 'A0137B0137 1kHz A weighted RMS vs pressing number side a')
% plot_scatter2(A0137B0137, '1kHz', 'b', 'PressingNumber', 'A_L', 'A_R', 'A0137B0137 1kHz A weighted RMS vs pressing number side b')
% plot_scatter2(A0000B0000, '1kHz', 'a', 'PressingNumber', 'A_L', 'A_R', 'A0000B0000 1kHz A weighted RMS vs pressing number side a')
% plot_scatter2(A0000B0000, '1kHz', 'b', 'PressingNumber', 'A_L', 'A_R', 'A0000B0000 1kHz A weighted RMS vs pressing number side b')


% plot_scatter2(A0137B0137, '1kHz', 'a', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'A0137B0137 1kHz CCIR weighted RMS vs pressing number side a')
% plot_scatter2(A0137B0137, '1kHz', 'b', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'A0137B0137 1kHz CCIR weighted RMS vs pressing number side b')
% plot_scatter2(A0000B0000, '1kHz', 'a', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'A0000B0000 1kHz CCIR weighted RMS vs pressing number side a')
% plot_scatter2(A0000B0000, '1kHz', 'b', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'A0000B0000 1kHz CCIR weighted RMS vs pressing number side b')


plot_stereohistogram(BOTH, 'quiet','a','RMS_L','RMS_R','Both pressings RMS levels in the quiet tracks side a','test')

plot_scatter2(BOTH, 'quiet', 'b', 'PressingNumber', 'RMS_R', 'RMS_L', 'test')
plot_scatter2(BOTH, '1kHz', 'b', 'PressingNumber', 'minMouldSteamIn_F', 'maxMouldSteamIn_F', 'test2')

%~~~~~~~~~~~~~ 1 kHz Ends ~~~~~~~~~~~~~%

% Tbl1 = A0000B0000(strcmp(Tbl.track,'quiet'),:);
% Tbl1 = Tbl1(strcmp(Tbl1.side,'a'),:);

% Tbl2 = A0137B0137(strcmp(A0137B0137.track,'quiet'),:);
% Tbl2 = Tbl2(strcmp(Tbl2.side,'a'),:);

% Tbl1.StamperHits = Tbl1.PressingNumber;
% Tbl2.StamperHits = Tbl2.PressingNumber + 197;

% x1 = Tbl1.StamperHits;
% x2 = Tbl2.StamperHits;

% stamperhits = [x1,x2];

% y1 = A0000B0000.RMS_L;
% y2 = A0137B0137.RMS_R;

% plot_scatterraw

function stats = get_stats(Tbl, trackname, side, column)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, column));
    X = table2array(Tbl(:,colx));

    stats = datastats(X);

end

function plot_barchart2(plotnum, Tbl, trackname, side, x1, x2,titlestring, varname, filename)
    cols = Tbl.Properties.VariableNames;
    pressruns = unique(Tbl.pressing);

    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx1 = find(ismember(cols, x1));
    colx2 = find(ismember(cols, x2));
    % X = table2array(Tbl(:,colx));

    X1 = [];
    X2 = [];
    for i = 1:length(pressruns)
        tbl = Tbl(strcmp(Tbl.pressing,pressruns{i}),:);
        x1 = mean(table2array(tbl(:,colx1)));
        x2 = mean(table2array(tbl(:,colx2)));
        X1 = [X1, x1];
        X2 = [X2, x2];
    end


    figure(plotnum); grid on; hold on;
    H = bar([X1, X2], 'LineWidth', 2);

    hold on;
    % H = bar(X1, 'LineWidth', 2)
    H(1).FaceColor = [0.6 0.6 0.6];
    H(2).FaceColor = [.9 .9 .9];
    set(gca,'xticklabel',pressruns)
    ax=gca;
    ax.FontSize=8;
    ax.XTick = (1:length(pressruns))   %THIS WAY, YOU SET HOW MANY XTICKS YOU WANT FOR YOUR XTICKLABELS
    xtickangle(45)
    ylabel(varname)
    title('RMS noise in quiet track')
    saveas(figure(plotnum),'RMSquiet.png')
end

function plot_barchart(plotnum, Tbl, trackname, side, x, titlestring, varname, filename)
    cols = Tbl.Properties.VariableNames;
    pressruns = unique(Tbl.pressing);

    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, x));
    % X = table2array(Tbl(:,colx));

    X = [];
    for i = 1:length(pressruns)
        tbl = Tbl(strcmp(Tbl.pressing,pressruns{i}),:);
        x = mean(table2array(tbl(:,colx)));
        X = [X, x];
    end


    figure(plotnum); grid on; hold on;
    H = bar(X, 'LineWidth', 2)
    H(1).FaceColor = [0.6 0.6 0.6];
    % H(2).FaceColor = [.9 .9 .9];
    set(gca,'xticklabel',pressruns)
    ax=gca;
    ax.FontSize=8;
    ax.XTick = (1:length(pressruns))   %THIS WAY, YOU SET HOW MANY XTICKS YOU WANT FOR YOUR XTICKLABELS
    xtickangle(45)
    ylabel(varname)
    title('RMS noise in quiet track')
    saveas(figure(plotnum),'RMSquiet.png')
end

function plot_histogram(plotnum, Tbl, trackname, side, x, titlestring, varname, binlims, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, x));
    X = table2array(Tbl(:,colx));

    figure(plotnum); grid on; hold on;
    histogram(X,50,'BinLimits',binlims)

    ylabel('number of records')
    varname
    xlabel(varname)
    title(titlestring)
    plotname = strcat('plots/', filename,'.png');
    saveas(figure(plotnum), plotname)
end


function plot_scatter(Tbl, trackname, side, x, y,titlestring, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, x));
    coly = find(ismember(cols, y));
    X = table2array(Tbl(:,colx));
    Y = table2array(Tbl(:,coly));

    % figure(plotnum);  
    fig = figure('Visible', 'off');
    scatter(X,Y,'ko');
    grid on; 
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end

function plot_scatter2(Tbl, trackname, side, x, y1, y2, titlestring, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, x));
    coly1 = find(ismember(cols, y1));
    coly2 = find(ismember(cols, y2));
    X = table2array(Tbl(:,colx));
    Y1 = table2array(Tbl(:,coly1));
    Y2 = table2array(Tbl(:,coly2));

    % figure(plotnum);  
    fig = figure('Visible', 'off');
    scatter(X,Y1,'ko');
    grid on; hold on;
    scatter(X,Y2,'kx');
    legend({'left channel', 'right channel'});
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end

function plot_scatterraw(X, Y, titlestring, filename)
    fig = figure('Visible', 'off');
    scatter(X,Y,'ko');
    grid on; 
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end

function plot_scatterraw2(X, Y, titlestring, filename)
    fig = figure('Visible', 'off');
    scatter(X,Y,'ko');
    grid on; 
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end



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
    
    % lower_binL = statsL.mean - 10;
    % lower_binR = statsR.mean - 10;
    % upper_binL = statsL.mean + 10;
    % upper_binR = statsR.mean + 10;


    fig = figure('Visible', 'off');
    % histogram(data_L, 50,'BinLimits',[lower_binL,upper_binL], 'facecolor',[0.6 0.6 0.6])
    histogram(data_L, 50, 'facecolor',[0.6 0.6 0.6])
    hold on; grid on;
    histogram(data_R,50, 'facecolor',[0.3 0.3 0.3])
    % histogram(data_R,50,'BinLimits',[lower_binL,upper_binL], 'facecolor',[0.3 0.3 0.3])
    title(titlestring)
    legend('left channel', 'right channel')
    dim = [0.2 0.5 0.3 0.3];
    str = {strcat('left mean :',num2str(statsL.mean)),strcat('left std :',num2str(statsL.std)),strcat('right mean :',num2str(statsR.mean)),strcat('right std :',num2str(statsR.std))};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor', 'white');
    plotname = strcat('plots/',filename,'.png');
    saveas(fig, plotname);
end