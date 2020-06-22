close all; %clear all; %clc;






% try 
% addpath('D:\Code\vinyl-research\matlab_code\audio_functions\')

% addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\')
% folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\';

% data_folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\'

% files = dir('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\*.wav')
% catch
    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/')
    % files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/*.wav')
% end
% % ~~~~ WINDOWS ~~~~ %
% addpath('D:\Code\vinyl-research\matlab_code\')
% addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
% addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\')
% addpath('E:\audio_files\A0000B0000\')


%~~~~~ TESTING ~~~~~%
% folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\';
%~~ TESTING ENDS ~~~%

% folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0137B0137\';
% addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
% folder = ('/Volumes/AUDIOBANK/audio_files/A0137B0137/')


% % ~~~~ WINDOWS ENDS ~~~~ %
% % ~~~~ MAC  ~~~~ %


if ismac() == true
    addpath('/Users/cz/Code/vinyl-research/matlab_code')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Wow')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
    addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')
    folder = ('/Volumes/AUDIOBANK/audio_files/A0000B0000/')


end 
if ispc() == true
    addpath('D:\Code\vinyl-research\matlab_code')
    addpath('D:\Code\vinyl-research\matlab_code\Common')
    addpath('D:\Code\vinyl-research\matlab_code\Wow')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin')    
    addpath('E:\audio_files\A0000B0000')
    folder = ('E:\audio_files\A0000B0000\')

end

% pressingID = 'A0137B0137';
pressingID = 'A0000B0000';


disp(['loading folder...:', folder])
% files = dir(strcat(folder,'*.wav'))
files = dir(fullfile(folder,'*.wav'))


% AudioTableHeaders = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};
AudioTableHeaders = {'pressing','record', 'side','track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicksa_L', 'commonclicksa_R','commonclicksb_L', 'commonclicksb_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};


% check if there is already a csv file to append to 
try
    disp('trying...')
    strcat(folder,strcat(pressingID,'-AudioTable.csv'))
    AudioTable = readtable(strcat(folder,strcat(pressingID,'-AudioTable.csv')))
catch
    disp('csv file not found, creating one...')
    AudioTable  = cell2table(cell(0,length(AudioTableHeaders)), 'VariableNames', AudioTableHeaders);
end
for i = (1:length(files)) %%loop through records
    % AudioTable
    % AudioTable.record
    % i
    filename = files(i).name;
    files(i);
    % filename = filename(3:end);
    disp(['opening file...:', filename])
    % filename(19:21)
    % AudioTable.record
    % filename(19:21)
    % num2str(AudioTable.record)
    % ismember(filename(19:21), num2str(AudioTable.record))
    % disp('ismember')
    % ismember(str2num(filename(19:21)), (AudioTable.record))
    % ismember(str2num(filename(22)), (AudioTable.side))
    
    if ismember(filename, AudioTable.record)
        disp('record already processed...')
        continue
    end
    
    file = strcat(files(i).folder,'/',filename);
    % file = strcat(files(i).folder,'\',filename)

    % pressid = 0; 
    date_recorded = 0;
    pressing = 0;
    top_stamper = 0;
    top_hits = 0;
    bottom_stamper = 0;
    bottom_hits = 0;
    record = 0;
    side  = 0;
    track  = 0;
    %     %STRIP RELEVANT INFO FROM NAME 
%     date_recorded = (filename(1:6));
%     % date_recorded = date_recorded{1};
%     record = str2num(filename(19:21));
%     % record = filename;
%     top_stamper = filename(8);
%     pressing = filename(8:16);
%     top_hits = str2num(filename(9:12)) + record;
%     bottom_stamper = filename(13);
%     bottom_hits = str2num(filename(14:17)) + record;
%     side = filename(end-4);

    % PressingNumber 
    % RecordID 
    % pressing 
    % RecordNumber 

    recordid = str2num(filename(19:21));

    % pressid = RecordNumbers.pressing(strcmp(RecordNumbers.RecordID,recordid),:)
    % pressing = RecordNumbers(RecordNumbers.RecordID == recordid, :);
    % pressing = pressing.pressing{1};
    pressing = 'A0000B0000';
    record = filename;
    side = filename(end-4);

    disp([strcat('...pressing:', pressing)])
    disp([strcat('...recordid:', recordid)])
    disp([strcat('...record:', record)])
    disp([strcat('...side:', side)])

    infoCell = {pressing, record, side};
    AudioOutput = RecordProcess(file);

    infoCell
    % cell2table([infoCell, AudioOutput], 'VariableNames', AudioTableHeaders)
    for j = (1:length(AudioOutput))
        AudioTable = [AudioTable; cell2table([infoCell, AudioOutput(j,:)], 'VariableNames', AudioTableHeaders)]
    end
    %append audio output to info cell array
    disp('SAVING CSV')
    writetable(AudioTable, strcat(folder,pressingID,'-AudioTable.csv'));
end





% files = dir(strcat(folder,'*.wav'));
% % % ~~~~ MAC ENDS ~~~~ %

% % addpath('E:\audio_files\A0000B0000\')
% % folder = ('E:\audio_files\A0000B0000\')

% disp(['loading folder...:', folder])

% files = dir(strcat(folder,'*.wav'));

% files

% % csvdata = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};
% csvdata = {'record', 'side','track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};


% pressingID = files.folder;
% pressingID = pressingID(end-9:end);

% % files.folder(1)
% try
%     T = readtable(strcat(folder,pressingID,'.csv'));
%     % T = table2cell(T);
% catch
%     disp('csv file not found, creating one...')
%     T  = cell2table(cell(0,length(csvdata)), 'VariableNames', csvdata);
%     % T = cell2table(empt,'VariableNames',csvdata(1,:))
%     % writetable(T,strcat(folder,pressingID,'.csv'));
% end


% for i = (1:length(files))
%     filename = files(i).name;
%     disp(['opening file...:', filename])

%     % T.record
%     if isempty(T.record)
%     elseif ismember(str2num(filename(19:21)), T.record)
%         disp('record already processed...')
%         continue
%     end
%     % disp('T{7}')
%     % T
%     % % if isempty(T{8,:})
%     % if ismember(str2num(filename(19:21)), [T{:,7}])
%     %     disp('record already processed...')
%     %     continue
%     % end

%     file = strcat(files(i).folder,'/',files(i).name);

%     date_recorded = 0;
%     pressing = 0;
%     top_stamper = 0;
%     top_hits = 0;
%     bottom_stamper = 0;
%     bottom_hits = 0;
%     record = 0;
%     side  = 0;
%     track  = 0;
    
%     %STRIP RELEVANT INFO FROM NAME 
%     date_recorded = (filename(1:6));
%     % date_recorded = date_recorded{1};
%     record = str2num(filename(19:21));
%     % record = filename;
%     top_stamper = filename(8);
%     pressing = filename(8:16);
%     top_hits = str2num(filename(9:12)) + record;
%     bottom_stamper = filename(13);
%     bottom_hits = str2num(filename(14:17)) + record;
%     side = filename(end-4);

%     top_hits = num2str(top_hits);
%     bottom_hits = num2str(bottom_hits);
%     record = num2str(record);

%     disp([strcat('...date_recorded:', date_recorded)])
%     disp([strcat('...pressing:', pressing)])
%     disp([strcat('...:record', record)])
%     disp([strcat('...top_stamper:', top_stamper)])
%     disp([strcat('...top_hits:', top_hits)])
%     disp([strcat('...bottom_stamper:', bottom_stamper)])
%     disp([strcat('...bottom_hits:', bottom_hits)])
%     disp([strcat('...side:', side)])


%     output = RecordProcess(file);
%     numrec = size(output);
    
%     for i=(1:numrec)
%         track = output(i,1);
%         lagdiff = output(i,2);
%         normalization_L = output(i,3);
%         normalization_R = output(i,4);
%         RMS_L = output(i,5);
%         RMS_R = output(i,6);
%         clicks_L = output(i,7);
%         clicks_R = output(i,8);
%         commonclicks_L = output(i,9);  
%         commonclicks_R = output(i,10);
%         THD_L = output(i,11);
%         THD_R = output(i,12);
%         wow_L = output(i,13);
%         wow_R = output(i,14);
%         stereo_bleed = output(i,15);
      
%         [track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, commonclicks_L, commonclicks_R  THD_L, THD_R, wow_L, wow_R, stereo_bleed];
        
%         % csv_towrite = [date_recorded, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, 
%         csv_towrite = [record, side ,track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, commonclicks_L, commonclicks_R, THD_L, THD_R, wow_L, wow_R, stereo_bleed];
%         T
%         T{:,end+1} = csv_towrite

%     end
%     disp('SAVING CSV')
%     % T = cell2table(csvdata(2:end,:),'VariableNames',csvdata(1,:));
%     writetable(T,strcat(folder,'A0137B0137.csv'));

% end

% % dlmwrite('test.csv',M,'delimiter',',');
% % N = randn(4,4);
% % dlmwrite('test.csv',N,'delimiter',',','-append');


