close all; %clear all; %clc;

% try 
addpath('D:\Code\vinyl-research\matlab_code\audio_functions\')

addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\')
folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\';

% files = dir('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\*.wav')
% catch
    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/')
    % files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/*.wav')
% end

% addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')
% folder = ('/Volumes/AUDIOBANK/audio_files/A0000B0000/')
% files = dir(strcat(folder,'*.wav'))

% addpath('E:\audio_files\A0000B0000\')
% folder = ('E:\audio_files\A0000B0000\')
disp(['loading folder...:', folder])

files = dir(strcat(folder,'*.wav'));

files

csvdata = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};



pressingID = files.folder;
pressingID = pressingID(end-9:end);

% files.folder(1)
try
    T = readtable(strcat(folder,pressingID,'.csv'));
    % T = table2cell(T);
catch
    disp('csv file not found, creating one...')
    T  = cell2table(cell(0,length(csvdata)), 'VariableNames', csvdata);
    % T = cell2table(empt,'VariableNames',csvdata(1,:))
    % writetable(T,strcat(folder,pressingID,'.csv'));
end


for i = (1:length(files))
    filename = files(i).name;
    disp(['opening file...:', filename])

    % T.record
    if isempty(T.record)
    elseif ismember(str2num(filename(19:21)), T.record)
        disp('record already processed...')
        continue
    end
    % disp('T{7}')
    % T
    % % if isempty(T{8,:})
    % if ismember(str2num(filename(19:21)), [T{:,7}])
    %     disp('record already processed...')
    %     continue
    % end

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
    % date_recorded = date_recorded{1};
    record = str2num(filename(19:21));

    top_stamper = filename(8);
    pressing = filename(8:16);
    top_hits = str2num(filename(9:12)) + record;
    bottom_stamper = filename(13);
    bottom_hits = str2num(filename(14:17)) + record;
    side = filename(22);

    top_hits = num2str(top_hits);
    bottom_hits = num2str(bottom_hits);
    record = num2str(record);

    disp([strcat('...date_recorded:', date_recorded)])
    disp([strcat('...pressing:', pressing)])
    disp([strcat('...:record', record)])
    disp([strcat('...top_stamper:', top_stamper)])
    disp([strcat('...top_hits:', top_hits)])
    disp([strcat('...bottom_stamper:', bottom_stamper)])
    disp([strcat('...bottom_hits:', bottom_hits)])
    disp([strcat('...side:', side)])

    % try 
        output = RecordProcess(file);
    % catch e
    %     disp(strcat('**CRASHED** recordProcess(''',file,''')'))
    %     fprintf(1,'The identifier was:\n%s',e.identifier);
    %     fprintf(1,'There was an error! The message was:\n%s',e.message);
    %     break
    % end
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
        commonclicks_L = output(i,9);  
        commonclicks_R = output(i,10);
        THD_L = output(i,11);
        THD_R = output(i,12);
        wow_L = output(i,13);
        wow_R = output(i,14);
        stereo_bleed = output(i,15);
      
        [track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, commonclicks_L, commonclicks_R  THD_L, THD_R, wow_L, wow_R, stereo_bleed];
        

        % lagdiff = lagdiff{1,1};
        % normalization_L = normalization_L{1,1};
        % normalization_R = normalization_R{1,1};
        % RMS_L = RMS_L{1,1};
        % RMS_R = RMS_R{1,1};
        % clicks_L = clicks_L{1,1};
        % clicks_R = clicks_R{1,1};
        % THD_L = THD_L{1,1};
        % THD_R = THD_R{1,1};
        % wow_L = wow_L{1,1};
        % wow_R = wow_R{1,1};
        % stereo_bleed = stereo_bleed{1,1};

        % disp([strcat('track...:', track{1})])
        % disp(['normalization_L...:', num2str(normalization_L)])
        % disp(['normalization_R...:', num2str(normalization_R)])
        % disp(['RMS_L...:', num2str(RMS_L)])
        % disp(['RMS_R...:', num2str(RMS_R)])
        % disp(['clicks_L...:', num2str(clicks_L)])
        % disp(['clicks_R...:', num2str(clicks_R)])
        % disp(['THD_L...:', num2str(THD_L)])
        % disp(['THD_R...:', num2str(THD_R)])
        % disp(['wow_L...:', num2str(wow_L)])
        % disp(['wow_R...:', num2str(wow_R)])
        % disp(['stereo_bleed...:', num2str(stereo_bleed)])
        % csvdata(end+1,:) = [date_recorded, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, record, side, track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow_L, wow_R, stereo_bleed];
        csv_towrite = [date_recorded, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, record, side, track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, commonclicks_L, commonclicks_R, THD_L, THD_R, wow_L, wow_R, stereo_bleed];
        
        % numbers = randi(9, 10, 1);
        % num_str = num2str(numbers);
        % num_cell = mat2cell(num_str, ones(10, 1), 1);

        % csvdata(1);
        % for k = (1:length(csvdata(end,:)))
        %     if iscell(csvdata(end,k))
        %         class(csvdata(end,k))
        %         tempcsv = csvdata(end,k)
        %         T.date_recorded(end)
        %         % csvdata(end,k) = tempcsv{1};
        %     end
        % end
        % out=cellfun(@num2str,csvdata, 'UniformOutput', false)
        % T
        % cell2table(csvdata(end,:))
        % T = [ T ; cell2table(csvdata(end,:),'VariableNames',csvdata(1,:)) ];
        T
        T{:,end+1} = csv_towrite
        % T_towrite = (cell2table(csv_towrite,'VariableNames', T.Properties.VariableNames))
        % date_recorded
        % class(date_recorded)
        % % date_recorded{1}
        % T_towrite.date_recorded{1} = str2num(date_recorded);
        % T_towrite.top_hits{1} = str2num(top_hits);
        % T_towrite.bottom_hits{1} = str2num(bottom_hits);
        % T_towrite.record{1} = str2num(record);

        % T_towrite.lagdiff{1} = str2double(lagdiff);
        % T_towrite.normalization_L{1} = str2double(normalization_L);
        % T_towrite.normalization_R{1} = str2double(normalization_R);
        % T_towrite.RMS_L{1} = str2double(RMS_L);
        % T_towrite.RMS_R{1} = str2double(RMS_R);
        % T_towrite.clicks_L{1} = str2double(clicks_L);
        % T_towrite.clicks_R{1} = str2double(clicks_R);
        % T_towrite.THD_L{1} = str2double(THD_L);
        % T_towrite.THD_R{1} = str2double(THD_R);
        % T_towrite.wow_L{1} = str2double(wow_L);
        % T_towrite.wow_R{1} = str2double(wow_R);
        % T_towrite.stereo_bleed{1} = str2double(stereo_bleed);



        % % T_towrite. = str2num();

        % % T_towrite. = str2num(T_towrite.);
        % T_towrite
        % T.date_recorded
        % T_towrite.date_recorded
        % T = [ T ; cell2table(csv_towrite,'VariableNames',T.Properties.VariableNames) ]
        % csv_towrite
        % csv_towrite = cell2mat(csv_towrite)
        % dlmwrite(strcat(folder, pressingID,'.csv'),csv_towrite.','-append');
        % % cell2csv('A0000B0000.csv',csv_towrite,',')
     end

end

% dlmwrite('test.csv',M,'delimiter',',');
% N = randn(4,4);
% dlmwrite('test.csv',N,'delimiter',',','-append');

% T = cell2table(csvdata(2:end,:),'VariableNames',csvdata(1,:))
% writetable(T,'myDataFile.csv');

function cell2csv(filename,cellArray,delimiter)
    % Writes cell array content into a *.csv file.
    % 
    % CELL2CSV(filename,cellArray,delimiter)
    %
    % filename      = Name of the file to save. [ i.e. 'text.csv' ]
    % cellarray    = Name of the Cell Array where the data is in
    % delimiter = seperating sign, normally:',' (it's default)
    %
    % by Sylvain Fiedler, KA, 2004
    % modified by Rob Kohr, Rutgers, 2005 - changed to english and fixed delimiter
    if nargin<3
        delimiter = ',';
    end
    
    datei = fopen(filename,'w');
    for z=1:size(cellArray,1)
        for s=1:size(cellArray,2)
            
            var = eval(['cellArray{z,s}']);
            
            if size(var,1) == 0
                var = '';
            end
            
            if isnumeric(var) == 1
                var = num2str(var);
            end
            
            fprintf(datei,var);
            
            if s ~= size(cellArray,2)
                fprintf(datei,[delimiter]);
            end
        end
        fprintf(datei,'\n');
    end
    fclose(datei);
end