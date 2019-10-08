
addpath('E:\audio_files\A0000B0000\')

audioFile = ('E:\audio_files\A0000B0000\A0000B0000.csv')
pressFile = ('E:\audio_files\A0000B0000\oct10A0000B0000.csv')

% opts = detectImportOptions(audioFile)
% opts = setvartype(opts, 'char')
% audioData = readtable(audioFile, opts)

audioData = readtable(audioFile)
pressData = readtable(pressFile)


% cellfun(@(c) num2str, audioData.date_recorded, 'UniformInput', false)
audioData.date_recorded = string(audioData.date_recorded);

%get all the variable names in both csv's and coming into one big cell array
headers = [audioData.Properties.VariableNames, pressData.Properties.VariableNames];
%create new table with all the needed variables 

%for loop appending rows to the new table
Tbl  = cell2table(cell(0,length(headers)), 'VariableNames', headers);

pres = 1;
for aud = 1:length(audioData.record)
    %matchup by RECORDID and record
    % num2str(audioData.record)
    if strcmp(pres,audioData.record(aud)) == 0
        for pres = 1:length(pressData.RECORDID)
            if strcmp(pressData.RECORDID(pres), audioData.record(aud))
                pres = pres;
                break
            end
        end
    end
    %now aud should be the row in audio and pres the corresponding row in press 

    for i = (1:length(audioData.Properties.VariableNames))
        aud,i
        audioData.Properties.VariableNames{i}
        audioData.Properties.VariableNames(i)
        x = audioData(aud,audioData.Properties.VariableNames{i})
        x = Tbl(aud,audioData.Properties.VariableNames(i))
        
        disp('printing')
        Tbl(aud,audioData.Properties.VariableNames(i)) = (audioData(aud,audioData.Properties.VariableNames(i)))
    end
    % disp('audioData for loop done')
    
    Tbl
    disp('pres for loop')
    for i = (1:length(pressData.Properties.VariableNames))
        i
        Tbl(aud,i) = pressData(aud, i);%i + length(audioData.Properties.VariableNames));
    end
    % disp('pressData for loop done')

    % Tbl

    % pres = strfind(pressData.RECORDID, num2str(audioData.record(aud)))
    % pressData.RECORDID(pres)


    %make all track data into one table entry called tracks

    %append rows and columns
    %convert all columns that can be to numbers in both
    %if #NaN
    %if a column doesnt exist in one table replace with '' or NaN 



end



% Tbl = join(audioData,pressData)

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

% figure(1); hold on; grid on;
% title('noise levels in "silent" tracks')
% xlabel('necord Number')
% ylabel('RMS level (dB)')

% % scatter(quiet.record, quiet.RMS_L)
% % scatter(quiet.record, quiet.RMS_R)
% % scatter(quiet2.record, quiet2.RMS_L)
% % scatter(quiet2.record, quiet2.RMS_R)
% scatter(transition.record, transition.RMS_L)
% scatter(transition.record, transition.RMS_R)




% figure(2); hold on; grid on;
% title('number of clicks in "silent" tracks')
% xlabel('necord number')
% ylabel('number of clicks')
% scatter(transition.record, transition.clicks_L )
% scatter(transition.record, transition.clicks_R )
% % scatter(quiet.record, quiet.clicks_L )
% % scatter(quiet.record, quiet.clicks_R )
% % scatter(quiet2.record, quiet2.clicks_L )
% % scatter(quiet2.record, quiet2.clicks_R )


% figure(3); hold on; grid on;
% title('number of clicks vs. pressure in "silent" tracks')
% xlabel('pressure (tons)')
% ylabel('number of clicks')
% scatter(transition.record, transition.clicks_L)
% scatter(transition.record, transition.clicks_R )


% s1kHz = Tbl(strcmp(Tbl.track,'1kHz'),:);
