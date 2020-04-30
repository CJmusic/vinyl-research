
% %~~~~~~~   TESTING   ~~~~~~~%
        
% % file = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0137B0137\maxoverlap3a.wav';
% % recordProcessTest(file)
% %~~~~~~~ TESTING ENDS ~~~~~~~%

% close all; clear all; clc;


% disp('~~~~~~~~~~~~TESTING RECORDPROCESS~~~~~~~~~~~~')

% % % ~~~~ WINDOWS ~~~~ %
% % addpath('D:\Code\vinyl-research\matlab_code\')
% % addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
% % addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\')
% % addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000')
% % addpath('E:\audio_files\A0000B0000\')


% % folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0137B0137\';
% % folder = 'E:\audio_files\A0000B0000\';
% % folder = ('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\macvsacer\')

% % % ~~~~ WINDOWS ENDS ~~~~ %


% if ismac() == true
%     addpath('/Users/cz/Code/vinyl-research/matlab_code')
%     addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
%     addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
%     addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
%     folder = ('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
%     RecordNumbers = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137_RecordNumbers.csv')

% end 
% if ispc() == true
%     addpath('D:\Code\vinyl-research\matlab_code')
%     addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
%     addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin')    
%     addpath('E:\audio_files\A0137B0137')
%     folder = ('E:\audio_files\A0137B0137\')
%     RecordNumbers = readtable('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0137B0137\A0137B0137_RecordNumbers.csv')

% end



% pressingID = 'A0137B0137';


% disp(['loading folder...:', folder])
% % files = dir(strcat(folder,'*.wav'))
% files = dir(fullfile(folder,'*.wav'))


% % AudioTableHeaders = {'date_recorded', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicks_L', 'commonclicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};
% AudioTableHeaders = {'pressing','record', 'side','track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'A_L', 'A_R', 'CCIR_L', 'CCIR_R','clicks_L', 'clicks_R', 'commonclicksa_L', 'commonclicksa_R','commonclicksb_L', 'commonclicksb_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};


% % check if there is already a csv file to append to 
% try
%     disp('trying...')
%     strcat(folder,strcat(pressingID,'-AudioTable.csv'))
%     AudioTable = readtable(strcat(folder,strcat(pressingID,'-AudioTable.csv')))
% catch
%     disp('csv file not found, creating one...')
%     AudioTable  = cell2table(cell(0,length(AudioTableHeaders)), 'VariableNames', AudioTableHeaders);
% end
% for i = (1:length(files)) %%loop through records
%     filename = files(i).name;
%     files(i);
%     disp(['opening file...:', filename])

    
%     if ismember(filename, AudioTable.record)
%         disp('record already processed...')
%         continue
%     end
    
%     file = strcat(files(i).folder,'/',filename);
%     date_recorded = 0;
%     pressing = 0;
%     top_stamper = 0;
%     top_hits = 0;
%     bottom_stamper = 0;
%     bottom_hits = 0;
%     record = 0;
%     side  = 0;
%     track  = 0;
    
%     % PressingNumber 
%     % RecordID 
%     % pressing 
%     % RecordNumber 

%     % recordid = str2num(filename(1:end - 5));
%     recordid = 0;

%     % pressid = RecordNumbers.pressing(strcmp(RecordNumbers.RecordID,recordid),:)
%     % pressing = RecordNumbers(RecordNumbers.RecordID == recordid, :);
%     % pressing = pressing.pressing{1};

%     record = filename;
%     pressing = 'lacquer'
%     % side = filename(end-4);

%     disp([strcat('...pressing:', pressing)])
%     disp([strcat('...recordid:', recordid)])
%     disp([strcat('...record:', record)])
%     disp([strcat('...side:', side)])

%     infoCell = {pressing, record, side};
%     AudioOutput = lacquerprocess(file);

%     infoCell
%     % cell2table([infoCell, AudioOutput], 'VariableNames', AudioTableHeaders)
%     for j = (1:length(AudioOutput))
%         AudioTable = [AudioTable; cell2table([infoCell, AudioOutput(j,:)], 'VariableNames', AudioTableHeaders)]
%     end
%     %append audio output to info cell array
%     disp('SAVING CSV')
%     writetable(AudioTable, strcat(folder,pressingID,'-AudioTable.csv'));
% end

file = 'd:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/lacquer_recordings/Dec 20 - Test A Part one .wav'
offset = 12.6;

lacquerprocess(file, offset)


function output = lacquerprocess(file, offset)
            % file(length(file)-4);
            % if ismac() == true
            %     % %~~~~ MAC ~~~~%
            %     disp('MAC')
            %     if file(length(file)-4) == 'a'
 
                 
 
            %     elseif file(length(file)-4) == 'b'
              
                
            %     end

         
            % end 
            % if ispc() == true
            %     %~~~~ WINDOWS ~~~~%
            %     if file(length(file)-4) == 'a'
                    
            %     end
            %     if file(length(file)-4) == 'b'
                   
            %     end

            % end
        %~~~~~~~~~~~~~~~~~~ Reference info ~~~~~~~~~~~~~~~~%
    
            timestamps_ref = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
            % this is how many seconds each signal is according to Chris Muths track listing
            lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
            signal_names = {'leadin',    % 1
                            '1kHz',      % 2
                            '10kHz',     % 3
                            '100Hz',     % 4
                            'sweep',     % 5
                            'quiet',     % 6
                            '3150Hz',    % 7 
                            '1kHzL',     % 8
                            'sweepL',    % 9
                            '1kHzR',     % 10
                            'sweepR',    % 11
                            '1kHzV',     % 12
                            'sweepV',    % 13
            };
            timestamps =       [[0, 61],    % 1. 1 kHz
                                [61,91],    % 2. 10 kHz
                                [91,121],   % 3. 100 Hz
                                [121,159],  % 4. sweep
                                [159,180],  % 5. quiet
                                [180,245],  % 6. 3150 Hz
                                [245,267],  % 7. 1 kHz left
                                [267, 302], % 8. sweep left
                                [302, 325], % 9. 1 kHz right
                                [325, 361], % 10. sweep right
                                [361, 383], % 11. 1 kHz vertical
                                [383, 418]] % 12. sweep vertical
                                % [938, 950]];               
                                %% dont forget lead in and leadout
            timestamps = timestamps + offset;
    
    
    
    
        %~~~~~~~~~~~~~~~~~    LOAD FILE    ~~~~~~~~~~~~~~~~%
    
            [data, fs] = audioread(file);
            output = {};
    
        %~~~~~~~~~~~~~~~~~~~~~ LINE UP ~~~~~~~~~~~~~~~~~~~~~~~~% 
    
            lagdiff = offset*fs;
            timedata = (0:length(data)-1)/fs  + lagdiff/fs;
    
           
    
        %~~~~~~~~~~~~~~~~~~~~ NORMALIZATION ~~~~~~~~~~~~~~~~~~~~%
            t = 1;
            sigtime = timedata(floor(timestamps(t,1)*fs):floor(timestamps(t,2)*fs));
            sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);
    
            % figure(t)
            % plot(sigtime,sig)
    
    
            % floor(timestamps(t,1)*fs) - lagdiff 
            % floor(timestamps(t,2)*fs) - lagdiff
    
            sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
            sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
            sigRMS=rms(sig);
            normalization=sqrt(2)*sigRMS*40/7; %digital value of peak level
            data(:,1)=data(:,1)/normalization(1);% now normalized to 40cm/s peak    
            data(:,2)=data(:,2)/normalization(2);% now normalized to 40cm/s peak 
            normalization_L = normalization(1);
            normalization_R = normalization(2);
    
            disp(strcat('normalization_L... ', num2str(size(normalization_L))))
            disp(strcat('normalization_R... ', num2str(size(normalization_R))))
            
            for t = 1:length(signal_names)
                % for each track we need: 
                %  - RMS level
                %  - Clicks
                %  - THD
                
                % for unique tracks we need 
                %  - normalization
                %  - stereo bleed (still really unsure about this test)
                %  - wow and flutter 
                track_name = signal_names{t};
                disp(strcat('track  ...',track_name))
    
    
                csig = [];
                CLICKS_R = [];
                CLICKS_L = [];
                RMS_L = [];
                RMS_R = [];
                THD_L = [];
                THD_R = [];
                
                if t == 1
                    sig = data(1 : floor(timestamps(1,1)*fs) - lagdiff,:);
                    sigtime = timedata(1 : floor(timestamps(1,1)*fs) - lagdiff);  
    
                    refT = ref(1 : floor(timestamps(1,1)*fs) - lagdiff,:);
                elseif t == length(signal_names)
                    sig = data(floor(timestamps(end,2)*fs) - lagdiff : length(data),:);
                    sigtime = timedata(floor(timestamps(end,2)*fs) - lagdiff : length(data));  
    
                    refT = ref(floor(timestamps(end,2)*fs) - lagdiff : length(ref),:);
                else
                    sig = data(floor(timestamps(t-1,1)*fs) - lagdiff : floor(timestamps(t-1,2)*fs) - lagdiff,:);
                    sigtime = timedata(floor(timestamps(t-1,1)*fs) - lagdiff :floor(timestamps(t-1,2)*fs) - lagdiff);  
    
                    % floor(timestamps(t-1,1)*fs) - lagdiff
                    % floor(timestamps(t-1,2)*fs) - lagdiff
    
                    refT = ref(floor(timestamps(t-1,1)*fs) - lagdiff : floor(timestamps(t-1,2)*fs) - lagdiff,:);
                end
                
                [csig(:,1), CLICKS_L] = ClickDetect(sig(:,1));
                [csig(:,2), CLICKS_R] = ClickDetect(sig(:,2));



               
                [~, REFSa_L] = ClickDetect(refT(:,1));
                [~, REFSa_R] = ClickDetect(refT(:,2));
                [~, REFSb_L] = ClickDetect(refT(:,1));
                [~, REFSb_R] = ClickDetect(refT(:,2));
    
                % need to do the reference here by track 
                commonclicksa_L = CommonClicks(CLICKS_L, REFSa_L);
                commonclicksa_R = CommonClicks(CLICKS_R, REFSa_R);
                commonclicksb_L = CommonClicks(CLICKS_L, REFSb_L);
                commonclicksb_R = CommonClicks(CLICKS_R, REFSb_R);
    
                
                clicks_L = length(CLICKS_L);
                clicks_R = length(CLICKS_R); 
    
                RMS_L = 20.0*log10(rms(csig(:,1)));
                RMS_R = 20.0*log10(rms(csig(:,2)));
    
                Aw = audio_Aweighting(csig(:,1));
                CCIRw = audio_CCIRweighting(csig(:,1));
    
                A_L = 20.0*log10(rms_response(Aw(1,:)));
                A_R = 20.0*log10(rms_response(Aw(2,:)));
    
                CCIR_L = 20.0*log10(avg_response(CCIRw(1,:)));
                CCIR_R = 20.0*log10(avg_response(CCIRw(2,:)));
    
                %***    DEBUG    ***%
                % figure(t); grid on;
                % plot(csig)
                % title(track_name)
                %*** DEBUG ENDS ***%
    
    
                THD_L = thd(csig(:,1),fs);
                THD_R = thd(csig(:,2),fs);
    
    
                if ismember(signal_names(t), {'1kHzL', 'sweepL', '1kHzL2', 'sweepL2'})             
                    stereo_bleed = StereoBleed(csig,1);
                elseif ismember(signal_names(t), {'1kHzR', 'sweepR', '1kHzR2', 'sweepR2'}) 
                    stereo_bleed = StereoBleed(csig,2);
                else
                    stereo_bleed = 0;
                end
    
                if ismember(signal_names(t), {'3150Hz', '3150Hz2'})
                    wow_L = WowFlutter(csig(:,1));
                    wow_R = WowFlutter(csig(:,2));
                else 
                    wow_L = 0;
                    wow_R = 0;
                end
              
                track = signal_names(t);
         
    
    
                output = [output; track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, A_L, A_R, CCIR_L, CCIR_R, clicks_L, clicks_R, commonclicksa_L, commonclicksa_R ,commonclicksb_L, commonclicksb_R, THD_L, THD_R, wow_L, wow_R, stereo_bleed];
            end
    end 