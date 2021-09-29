clear all; close all;
Tbl = readtable( '/Users/cz/Code/vinyl-research/matlab_code/A0137B0137/A0137B0137.csv');


% plot_scatter2(Tbl, 'quiet', 'a', 'PressingNumber', 'RMS_L', 'RMS_R', 'Pressing Number vs Unweighted RMS side a')
% plot_scatter2(Tbl, 'quiet', 'a', 'PressingNumber', 'A_L', 'A_R', 'Pressing Number vs A weighted RMS side a')
% plot_scatter2(Tbl, 'quiet', 'a', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'Pressing Number vs CCIR weighted RMS side a')

% plot_scatter2(Tbl, 'quiet', 'b', 'PressingNumber', 'RMS_L', 'RMS_R', 'Pressing Number vs Unweighted RMS side b')
% plot_scatter2(Tbl, 'quiet', 'b', 'PressingNumber', 'A_L', 'A_R', 'Pressing Number vs A weighted RMS side b')
% plot_scatter2(Tbl, 'quiet', 'b', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'Pressing Number vs CCIR weighted RMS side b')

% plot_scatter2(Tbl, 'quiet', 'a', 'PressingNumber', 'clicks_L', 'clicks_R', 'Pressing Number vs Unweighted RMS side a')
% plot_scatter2(Tbl, 'quiet', 'b', 'PressingNumber', 'CCIR_L', 'CCIR_R', 'Pressing Number vs CCIR weighted RMS side b')


% figure(4)
% plot_scatter2special(Tbl, 'quiet', 'a', 'PressingNumber', 'clicks_L', 'clicks_R', 'Pressing Number vs Clicks side a')
% figure(5)
% plot_scatter2(Tbl, 'quiet', 'b', 'PressingNumber', 'clicks_L', 'clicks_R', 'Pressing Number vs Clicks side b')

% figure(10)
% plot_histogram(10, Tbl, 'quiet', 'a', 'A_L', 'clicks histogram left channel', 'clicks', [-58 -52])



% plot_scatter6(Tbl, 'quiet', 'quiet2', 'a',  'PressingNumber', 'RMS_L', 'RMS_R', 'Unweighted RMS levels per record side a', 'A0000B0000quiettracksRMSa')
% plot_scatter6(Tbl, 'quiet', 'quiet2', 'a',  'PressingNumber', 'A_L', 'A_R', 'A weighted RMS levels per record side a', 'A0000B000quiettracksAa')
% plot_scatter6(Tbl, 'quiet', 'quiet2', 'a',  'PressingNumber', 'CCIR_L', 'CCIR_R', 'CCIR weighted RMS levels per record side a', 'A0000B000quiettracksCCIRa')

% plot_scatter6(Tbl, 'quiet', 'quiet2', 'b',  'PressingNumber', 'RMS_L', 'RMS_R', 'Unweighted RMS levels per record side b', 'A0000B0000quiettracksRMSb')
% plot_scatter6(Tbl, 'quiet', 'quiet2', 'b',  'PressingNumber', 'A_L', 'A_R', 'A weighted RMS levels per record side b', 'A0000B000quiettracksAb')
% plot_scatter6(Tbl, 'quiet', 'quiet2', 'b',  'PressingNumber', 'CCIR_L', 'CCIR_R', 'CCIR weighted RMS levels per record side b', 'A0000B000quiettracksCCIRb')


plot_clicks(Tbl, 'a','Total number of clicks in the first pressing side a')
plot_clicks(Tbl, 'b','Total number of clicks in the first pressing side b')


function plot_clicks(Tbl,side, titlestring)
    clicktracks = {{'100Hz'     }
    {'100Hz2'    }
    {'10kHz'     }
    {'10kHz2'    }
    {'1kHz'      }
    {'1kHz2'     }
    {'1kHzL'     }
    {'1kHzL2'    }
    {'1kHzR'     }
    {'1kHzR2'    }
    {'1kHzV'     }
    {'1kHzV2'    }
    {'3150Hz'    }
    {'3150Hz2'   }
    {'quiet'     }
    {'quiet2'    }
    {'sweep'     }
    {'sweep2'    }
    {'sweepL'    }
    {'sweepL2'   }
    {'sweepR'    }
    {'sweepR2'   }
    {'sweepV'    }
    {'sweepV2'   }
    {'transition'}};
    Tbl = Tbl(strcmp(Tbl.side,side),:);

    cols = Tbl.Properties.VariableNames;
    % clicks_L = cell2table(cell(0,length(cols)));
    % clicks_R = cell2table(cell(0,length(cols)));
    clicks_L = [];
    clicks_R = [];
    % for i = (1:length(clicktracks))
    Tbl(strcmp(Tbl.track,'leadout'),:) = [];
    Tbl(strcmp(Tbl.track,'leadin'),:) = [];
    PressingNumbers = unique(Tbl.PressingNumber);
    for i = (1:length(PressingNumbers))
        tbl = Tbl((Tbl.PressingNumber==PressingNumbers(i)),:);
        clicks_l = sum(tbl.clicks_L);
        clicks_r = sum(tbl.clicks_R);
        % for j = 1:length(clicktracks)
        %     tbl = Tbl(strcmp(Tbl.track,clicktracks{j}),:);
        %     clicks_l = tbl.clicks_L
        %     clicks_r = tbl.clicks_R
        % end
    
        clicks_L = [clicks_L; (clicks_l)];
        clicks_R = [clicks_R; (clicks_r)];
    
    end 
    pressing_number = PressingNumbers;
    size(pressing_number)
    size(clicks_L)
    size(clicks_R)
    % pressing_number = (1:137);
    fig = figure('Visible', 'off');
    scatter(pressing_number,clicks_L,'ko');
    grid on; hold on;
    scatter(pressing_number,clicks_R,'kx');
    legend({'left channel', 'right channel'});
    title(titlestring);
    % xlim([0,100])
    ylim([1,160])
    xlabel('record number')
    ylabel('number of clicks')
    plotname = strcat('plots/', titlestring,'.png');
    saveas(fig, plotname);

end

%~~~~~~~~~~~~~ Plot 
function plot_scatter6(Tbl, trackname1, trackname2, side, x, y1, y2, titlestring, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    Tbl1 = Tbl(strcmp(Tbl.track,trackname1),:);
    Tbl2 = Tbl(strcmp(Tbl.track,trackname2),:);
    % Tbl3 = Tbl(strcmp(Tbl.track,trackname3),:);
    colx = find(ismember(cols, x));
    coly1 = find(ismember(cols, y1));
    coly2 = find(ismember(cols, y2));
    X1 = table2array(Tbl1(:,colx));
    Y1L = table2array(Tbl1(:,coly1));
    Y1R = table2array(Tbl1(:,coly2));
    X2 = table2array(Tbl2(:,colx));
    Y2L = table2array(Tbl2(:,coly1));
    Y2R = table2array(Tbl2(:,coly2));
    % X3 = table2array(Tbl3(:,colx));
    % Y3L = table2array(Tbl3(:,coly1)); 
    % Y3R = table2array(Tbl3(:,coly2));

    % figure(plotnum);  
    fig = figure('Visible', 'off');
    scatter(X1,Y1L,'ko');
    grid on; hold on;
    scatter(X1,Y1R,'kx');
    scatter(X2,Y2L,'bo');
    grid on; hold on;
    scatter(X2,Y2R,'bx');
    % scatter(X3,Y3L,'go');
    % grid on; hold on;
    % scatter(X3,Y3R,'gx');
    legend({'quiet left', 'quiet right', 'quiet2 left', 'quiet2 right'});
    ylim([-58,-38])
    xlim([0,100])
    xlabel('record number')
    ylabel('Level [dB]')
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);
end



function plot_barchart2(plotnum, Tbl, trackname, side, x1, x2,titlestring, varname)
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




function plot_barchart(plotnum, Tbl, trackname, side, x, titlestring, varname)
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

function plot_histogram(plotnum, Tbl, trackname, side, x, titlestring, varname, binlims)
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
    plotname = strcat('plots/', titlestring,'.png');
    saveas(figure(plotnum), plotname)
end


function plot_scatter(Tbl, trackname, side, x, y,titlestring)
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
    plotname = strcat('plots/', titlestring,'.png');
    saveas(fig, plotname);

end

function plot_scatter2(Tbl, trackname, side, x, y1, y2, titlestring)
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
    plotname = strcat('plots/', titlestring,'.png');
    saveas(fig, plotname);

end

function plot_scatter2special(Tbl, trackname, side, x, y1, y2, titlestring)
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
    ylim([0,10])
    legend({'left channel', 'right channel'});
    title(titlestring);
    plotname = strcat('plots/', titlestring,'.png');
    saveas(fig, plotname);

end
