clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)



if ismac() == true
    addpath('/Users/cz/Code/vinyl-research/matlab_code')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
    addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Wow')
    folder = ('/Volumes/AUDIOBANK/audio_files/styluswear/')
end 
if ispc() == true
    addpath('D:\Code\vinyl-research\matlab_code')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin')    
    addpath('e:/audio_files/styluswear/')
    addpath('D:\Code\vinyl-research\matlab_code\A0137B0137')
    addpath('D:\Code\vinyl-research\matlab_code\Common\')
    folder = ('e:/audio_files/styluswear/')
end



pressingID = 'StylusWear';


disp(['loading folder...:', folder])
files = dir(fullfile(folder,'*.wav'))
AudioTableHeaders = {'pressing','record', 'side','track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicksa1_L', 'commonclicksa1_R','commonclicksb1_L', 'commonclicksb1_R', 'commonclicksa2_L', 'commonclicksa2_R','commonclicksb2_L', 'commonclicksb2_R',  'RMSclicks_L', 'RMSclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'centreholeoffset','stereo_bleed'};


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
    filename = files(i).name;
    files(i);
    disp(['opening file...:', filename])    
    if ismember(filename, AudioTable.record)
        disp('record already processed...')
        continue
    end
    
    file = strcat(files(i).folder,'/',filename);
    date_recorded = 0;
    pressing = 0;
    top_stamper = 0;
    top_hits = 0;
    bottom_stamper = 0;
    bottom_hits = 0;
    record = 0;
    side  = 0;
    track  = 0;
    recordid = str2num(filename(1:end - 5));

    pressing = 'styluswear';

    record = filename;
    side = filename(end-4);

    disp([strcat('...pressing:', pressing)])
    disp([strcat('...recordid:', recordid)])
    disp([strcat('...record:', record)])
    disp([strcat('...side:', side)])

    infoCell = {pressing, record, side};
    AudioOutput = RecordProcess(file);

    infoCell
    len = size(AudioOutput);
    len = len(1);
    for j = (1:len)
        AudioTable = [AudioTable; cell2table([infoCell, AudioOutput(j,:)], 'VariableNames', AudioTableHeaders)]
    end
    disp('SAVING CSV')
    writetable(AudioTable, strcat(folder,pressingID,'-AudioTable.csv'));
end



