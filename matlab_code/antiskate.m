clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)





H = bar([Tbl.AvgRMS_L(strcmp(AudioStats.track,'transition'),:), Tbl.AvgRMS_R(strcmp(AudioStats.track,'transition'),:)], 'LineWidth', 2)

H(1).FaceColor = [0.6 0.6 0.6];
H(2).FaceColor = [.9 .9 .9];
