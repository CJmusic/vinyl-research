
addpath('E:\audio_files\A0000B0000\')
dataFile = ('E:\audio_files\A0000B0000\')

dataFile = ('E:\audio_files\A0000B0000\A0000B0000.csv')


Tbl = readtable(dataFile);
% Tbl.track
Tbl.Properties.VariableNames
Tbl.track


transition = Tbl(strcmp(Tbl.track,'transition'),:);
quiet = Tbl(strcmp(Tbl.track,'transition'),:);
quiet2 = Tbl(strcmp(Tbl.track,'transition'),:);

figure(1); hold on; grid on;
scatter(quiet.record, quiet.RMS_L)
scatter(quiet.record, quiet.RMS_R)
scatter(quiet2.record, quiet2.RMS_L)
scatter(quiet2.record, quiet2.RMS_R)
scatter(transition.record, transition.RMS_L)
scatter(transition.record, transition.RMS_R)

title('RMS noise per record number')
  


figure(2); hold on; grid on;
scatter(transition.record, transition.clicks_L )
scatter(transition.record, transition.clicks_R )
scatter(quiet.record, quiet.clicks_L )
scatter(quiet.record, quiet.clicks_R )
scatter(quiet2.record, quiet2.clicks_L )
scatter(quiet2.record, quiet2.clicks_R )
title('clicks per record')


s1kHz = Tbl(strcmp(Tbl.track,'1kHz'),:);
