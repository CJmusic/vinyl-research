% audio_recordclass
% christopher zaworski
% last edit : march 30 2019 
%
% This is the matlab implementation fo the record class. 
%

classdef audio_recordclass < handle %inheriting handle allows methods to update properties
    properties 
        data = [];
        dataL = [];
        dataR = [];
        timeL = [];
        timeR = [];
        fs = 0;

        clicks = [];
        clicksL = [];
        clicksR = [];

        click_matrix = [];

        signalsL = []; %
        signalsR = []; %
        signals = {};
        
        %%% RECORD INFO 
        timestamps = [0, 2, 62, 92, 124, 160, 182, 248, 268, 306, 326, 364, 384, 419.5];
        % this is how many seconds each signal is according to Chris Muth's track listing
        lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
        signal_names = {'needledrop','leadin','1kHz', '10kHz', '100Hz', 'freqsweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition','leadin2','1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2', '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout2'};

        offset = 4.25; % as measured on /020818_A0000B0000/02072019_A0000B000r25-A.wav
        %time = (0:length(dataL)-1)/rec.fs;
        transition = 518.25; % as measured on /020818_A0000B0000/02072019_A0000B000r25-A.wav
        lagdiff = 0;
        timediff = 0; % this is a time diff calculated as compared to a reference 
        tracks;
    end % properties
    methods 
        function rec = audio_recordclass(file_path);
            disp('record class constructor')
            file_path
            [data, rec.fs] = audioread(file_path);
            rec.data = data;
            rec.dataL = data(:,1);
            rec.dataR = data(:,2);
            %rec.clicksL = audio_clickdetect(rec.dataL, rec.fs);
            %rec.clicksR = audio_clickdetect(rec.dataR, rec.fs);
        end

        function lagcorrect(rec, lagDiff);
            % this function should correct the record for any lagDiff calculated by another method: 
            % the method currently implemented is okay, I want to try circshift
            % also see the matlab function lag !!!! 
           time = (0:length(rec.data)); 
           if lagDiff > 0;
               % positive lagDiffL means that the data array is delayed compared to the reference
               cdata = rec.data(lagDiff + 1:end,:); 
               ctime = time - lagDiff/rec.fs; 
           elseif lagDiff < 0; 
               % negative lagDiffL means that the data array is ahead of the reference  
               cdataL = rec.data(1: end + lagDiff,:);   
               ctimeL = time - lagDiff/rec.fs;
           elseif lagDiff == 0; %then nothing needs to be corrected  
               cdata = rec.data; %no lag 
               ctime = time - lagDiff/rec.fs; 
           end

           
           rec.data = cdata;

           rec.time = ctime;
        end % function lagcorrect

        function process_tracks(rec);
            disp('Processing tracks')
            rec.offset;
            %TRACK LISTINGS ON TEST RECORDS 
            %00:02 : 1KHz@7cm-s lateral
            %01:02 : 10kHz @ -20dB
            %01:32 : 100Hz
            %02:04 : frequency sweep 20Hz-16kHz 
            %02:40 : quiet groove
            %03:02 : 3150 wow & flutter
            %04:08 : 1kHz left
            %04:28 : frequency sweep left
            %05:06 : 1kHz right
            %05:26 : frequency sweep right
            %06:04 : 1kHz vertical
            %06:24 : frequency sweep vertical

            rec.tracks = 0; % clear any previous tracks info
            rec.signals = {};

            timestamps  = rec.timestamps + rec.offset + rec.timediff;
            timestamps2 = timestamps + rec.transition + rec.offset + rec.timediff;
            
            % get the needledrop
            signalsL = rec.dataL(1:floor((timestamps(1))*rec.fs));
            signalsR = rec.dataR(1:floor((timestamps(1))*rec.fs));
            rec.signals{end + 1} = [signalsL, signalsR];
            
            % first set of signals on the disk
            for i = (1:length(timestamps)-1);
                signalsL = rec.dataL(floor(timestamps(i)*rec.fs):floor(timestamps(i+1)*rec.fs));
                signalsR = rec.dataR(floor(timestamps(i)*rec.fs):floor(timestamps(i+1)*rec.fs)); 
                rec.signals{end + 1} = [signalsL, signalsR];
            end
            
            % the extended silence section 'transition' between the two sets of signals 
            signalsL = rec.dataL(floor(timestamps(end)*rec.fs): ...
                                           floor((timestamps2(1))*rec.fs));
            
            signalsR = rec.dataR(floor(timestamps(end)*rec.fs): ...
                                           floor((timestamps2(1))*rec.fs));
            rec.signals{end + 1} = [signalsL, signalsR];
            
            % second set of signals on the disk
            for i = (1:length(timestamps2)-1);
                signalsL =rec.dataL(floor(timestamps2(i)*rec.fs):floor(timestamps2(i+1)*rec.fs));
                signalsR = rec.dataR(floor(timestamps2(i)*rec.fs):floor(timestamps2(i+1)*rec.fs));
                rec.signals{end + 1} = [signalsL, signalsR];
            end
            
            % get the rest of the file
            signalsL =rec.dataL(floor(timestamps2(end)*rec.fs):end);
            signalsR = rec.dataR(floor(timestamps2(end)*rec.fs):end);
            rec.signals{end + 1} = [signalsL, signalsR];
            
            rec.tracks = containers.Map(rec.signal_names, rec.signals);
            disp('printing tracks')
            rec.tracks
            disp('done processing tracks')

        end % function signals
        
        function clickdetect(rec);
            rec.clicks = audio_clickdetect(rec.data, rec.fs); 
        end % function detect_clicks

        function clicklineup(rec, clicks_ref)
            [rec.click_matrix, rec.lagdiff] = audio_clickmatrix(rec.clicks, clicks_ref);
            disp('lagdiff')
            size(rec.click_matrix)
            size(rec.lagdiff)
            rec.lagdiff
            rec.timediff = rec.lagdiff/rec.fs;
            rec.data = circshift(rec.data, rec.lagdiff);
            
        end % function click_lineup
    end % methods
end % recordclass
