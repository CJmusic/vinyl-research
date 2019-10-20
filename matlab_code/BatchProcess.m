close all; %clear all; %clc;

folder = '';
% try 
    % addpath('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\')
    % files = dir('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\click_testing\*.wav')
% catch
    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/')
    % files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/*.wav')
% end

% addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')
% folder = ('/Volumes/AUDIOBANK/audio_files/A0000B0000/')
% files = dir(strcat(folder,'*.wav'))

addpath('E:\audio_files\A0000B0000\')
folder = ('E:\audio_files\A0000B0000\')
disp(['loading folder...:', folder])

files = dir(strcat(folder,'*.wav'));

csvdata = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'clicks_L', 'clicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};



pressingID = files.folder;
pressingID = pressingID(end-9:end);

% files.folder(1)
try

    T = readtable(strcat(folder,pressingID,'.csv'));
    T = readtable(strcat(folder, pressingID,'.csv'));
catch
    disp('csv file not found, creating one...')
    T  = cell2table(cell(0,length(csvdata)), 'VariableNames', csvdata);
    % T = cell2table(empt,'VariableNames',csvdata(1,:))
    % writetable(T,strcat(folder,pressingID,'.csv'));
end


for i = (1:length(files))
    filename = files(i).name;
    disp(['opening file...:', filename])

    if isempty(T.record)
    elseif ismember(str2num(filename(19:21)), T.record)
        disp('record already processed...')
        continue
    end

    file = strcat(files(i).folder,'/',files(i).name);

    date_recorded = 0;
    pressing = 0;
    top_stamper = 0;
    top_hits = 0;
    bottom_stamper = 0;
    bottom_hits = 0;
    record = 0;
    side  = 0;
    track  = 0;
    
    %STRIP RELEVANT INFO FROM NAME 
    date_recorded = (filename(1:6));
    record = str2num(filename(19:21));

    top_stamper = filename(8);
    pressing = filename(8:16);
    top_hits = str2num(filename(9:12)) + record;
    bottom_stamper = filename(13);
    bottom_hits = str2num(filename(14:17)) + record;
    side = filename(22);

    disp([strcat('...date_recorded:', date_recorded)])
    disp([strcat('...pressing:', pressing)])
    disp([strcat('...:record', record)])
    disp([strcat('...top_stamper:', top_stamper)])
    disp([strcat('...top_hits:', top_hits)])
    disp([strcat('...bottom_stamper:', bottom_stamper)])
    disp([strcat('...bottom_hits:', bottom_hits)])
    disp([strcat('...side:', side)])

    try 
        output = recordProcess(file);
    catch
        disp(strcat('**CRASHED** recordProcess(''',file,''')'))
        break
    end
        numrec = size(output);
    for i=(1:numrec)
        track = output(i,1);
        lagdiff = output(i,2);
        normalization_L = output(i,3);
        normalization_R = output(i,4);
        RMS_L = output(i,5);
        RMS_R = output(i,6);
        clicks_L = output(i,7);
        clicks_R = output(i,8);
        THD_L = output(i,9);
        THD_R = output(i,10);
        wow_L = output(i,11);
        wow_R = output(i,12);
        stereo_bleed = output(i,13);
      
 
        

        lagdiff = lagdiff{1,1};
        normalization_L = normalization_L{1,1};
        normalization_R = normalization_R{1,1};
        RMS_L = RMS_L{1,1};
        RMS_R = RMS_R{1,1};
        clicks_L = clicks_L{1,1};
        clicks_R = clicks_R{1,1};
        THD_L = THD_L{1,1};
        THD_R = THD_R{1,1};
        wow_L = wow_L{1,1};
        wow_R = wow_R{1,1};
        stereo_bleed = stereo_bleed{1,1};

        disp([strcat('track...:', track{1})])
        disp(['normalization_L...:', num2str(normalization_L)])
        disp(['normalization_R...:', num2str(normalization_R)])
        disp(['RMS_L...:', num2str(RMS_L)])
        disp(['RMS_R...:', num2str(RMS_R)])
        disp(['clicks_L...:', num2str(clicks_L)])
        disp(['clicks_R...:', num2str(clicks_R)])
        disp(['THD_L...:', num2str(THD_L)])
        disp(['THD_R...:', num2str(THD_R)])
        disp(['wow_L...:', num2str(wow_L)])
        disp(['wow_R...:', num2str(wow_R)])
        disp(['stereo_bleed...:', num2str(stereo_bleed)])
        csvdata(end+1,:) = {date_recorded, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, record, side, track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow_L, wow_R, stereo_bleed};

        numbers = randi(9, 10, 1);
        num_str = num2str(numbers);
        num_cell = mat2cell(num_str, ones(10, 1), 1);

        % T = [ T ; cell2table(csvdata(end,:),'VariableNames',csvdata(1,:)) ];
        dlmwrite(strcat(folder, pressingID,'.csv'),cell2table(csvdata(end,:),'delimiter',',','-append'));
        % writetable(T,strcat(folder, pressingID,'.csv'));
    end

end

% dlmwrite('test.csv',M,'delimiter',',');
% N = randn(4,4);
% dlmwrite('test.csv',N,'delimiter',',','-append');

% T = cell2table(csvdata(2:end,:),'VariableNames',csvdata(1,:))
% writetable(T,'myDataFile.csv');