
addpath('E:\audio_files\A0000B0000\')

audioFile = ('E:\audio_files\A0000B0000\A0000B0000.csv')
pressFile = ('E:\audio_files\A0000B0000\oct10A0000B0000.csv')


audioData = readtable(audioFile);
pressData = readtable(pressFile)
% Tbl.track
% Tbl.Properties.VariableNames
% Tbl.track

Tbl = join(audioData,pressData)

% leadin = Tbl(strcmp(Tbl.track,'leadin'),:);
% s1kHz= Tbl(strcmp(Tbl.track,'1kHz'),:);
% s10kHz= Tbl(strcmp(Tbl.track,'10kHz'),:);
% s100Hz= Tbl(strcmp(Tbl.track,'100Hz'),:);
% sweep= Tbl(strcmp(Tbl.track,'sweep'),:);
% quiet= Tbl(strcmp(Tbl.track,'quiet'),:);
% s3150Hz= Tbl(strcmp(Tbl.track,'3150Hz'),:);
% s1kHzL= Tbl(strcmp(Tbl.track,'1kHzL'),:);
% sweepL= Tbl(strcmp(Tbl.track,'sweepL'),:);
% s1kHzR= Tbl(strcmp(Tbl.track,'1kHzR'),:);
% sweepR= Tbl(strcmp(Tbl.track,'sweepR'),:);
% s1kHzV= Tbl(strcmp(Tbl.track,'1kHzV'),:);
% sweepV= Tbl(strcmp(Tbl.track,'sweepV'),:);
% transition= Tbl(strcmp(Tbl.track,'transition'),:);
% s1kHz2= Tbl(strcmp(Tbl.track,'1kHz2'),:);
% s10kHz2= Tbl(strcmp(Tbl.track,'10kHz2'),:);
% s100Hz2= Tbl(strcmp(Tbl.track,'100Hz2'),:);
% sfreqsweep2= Tbl(strcmp(Tbl.track,'freqsweep2'),:);
% quiet2= Tbl(strcmp(Tbl.track,'quiet2'),:);
% s3150Hz2= Tbl(strcmp(Tbl.track,'3150Hz2'),:);
% s1kHzL2= Tbl(strcmp(Tbl.track,'1kHzL2'),:);
% sweepL2= Tbl(strcmp(Tbl.track,'sweepL2'),:);
% s1kHzR2= Tbl(strcmp(Tbl.track,'1kHzR2'),:);
% sweepR2= Tbl(strcmp(Tbl.track,'sweepR2'),:);
% s1kHzV2= Tbl(strcmp(Tbl.track,'1kHzV2'),:);
% sweepV2= Tbl(strcmp(Tbl.track,'sweepV2'),:);
% leadout= Tbl(strcmp(Tbl.track,'leadout'),:);

figure(1); hold on; grid on;
title('noise levels in "silent" tracks')
xlabel('necord Number')
ylabel('RMS level (dB)')

% scatter(quiet.record, quiet.RMS_L)
% scatter(quiet.record, quiet.RMS_R)
% scatter(quiet2.record, quiet2.RMS_L)
% scatter(quiet2.record, quiet2.RMS_R)
scatter(transition.record, transition.RMS_L)
scatter(transition.record, transition.RMS_R)




figure(2); hold on; grid on;
title('number of clicks in "silent" tracks')
xlabel('necord number')
ylabel('number of clicks')
scatter(transition.record, transition.clicks_L )
scatter(transition.record, transition.clicks_R )
% scatter(quiet.record, quiet.clicks_L )
% scatter(quiet.record, quiet.clicks_R )
% scatter(quiet2.record, quiet2.clicks_L )
% scatter(quiet2.record, quiet2.clicks_R )


figure(3); hold on; grid on;
title('number of clicks vs. pressure in "silent" tracks')
xlabel('pressure (tons)')
ylabel('number of clicks')
scatter(transition.record, transition.clicks_L)
scatter(transition.record, transition.clicks_R )


s1kHz = Tbl(strcmp(Tbl.track,'1kHz'),:);
