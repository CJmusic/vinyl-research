% audio_recordclass
% christopher zaworski
% last edit : march 30 2019 
%
% This is the matlab implementation fo the record class. 
%



classdef audio_recordclass < handle %inheriting handle allows methods to update properties
    properties 
        dataL = [];
        dataR = [];
        timeL = [];
        timeR = [];
        fs = 0;

        clicksL = [];
        clicksR = [];

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
        tracks;
    end % properties
    methods 
        function rec = audio_recordclass(file_path);
            disp('record class constructor')
            file_path
            [data, rec.fs] = audioread(file_path);
            rec.dataL = data(:,1);
            rec.dataR = data(:,2);
            %rec.clicksL = audio_clickdetect(rec.dataL, rec.fs);
            %rec.clicksR = audio_clickdetect(rec.dataR, rec.fs);
        end

        function lagcorrect(rec, lagDiffL, lagDiffR);

           timeL = (0:length(rec.dataL)); 
           if lagDiffL > 0;
               % positive lagDiffL means that the data array is delayed compared to the reference
               cdataL = rec.dataL(lagDiffL + 1:end,:); 
               ctimeL = timeL - lagDiffL/rec.fs; 
           elseif lagDiffL < 0; 
               % negative lagDiffL means that the data array is ahead of the reference  
               cdataL = rec.dataL(1: end + lagDiffL,:);   
               ctimeL = timeL - lagDiffL/rec.fs;
           elseif lagDiffL == 0; %then nothing needs to be corrected  
               cdataL = rec.dataL; %no lag 
               ctimeL = timeL - lagDiffL/rec.fs; 
           end


           timeR = (0:length(rec.dataR)); 
           if lagDiffR > 0;
               % positive lagDiffR means that the data array is delayed compared to the reference
               cdataR = rec.dataR(lagDiffR + 1:end,:); 
               ctimeR = timeR - lagDiffR/rec.fs; 
           elseif lagDiffR < 0; 
               % negative lagDiffR means that the data array is ahead of the reference  
               cdataR = rec.dataR(1: end + lagDiffR,:);   
               ctimeR = timeR - lagDiffR/rec.fs;
           elseif lagDiffR == 0; %then nothing needs to be corrected  
               cdataR = rec.dataR; %no lag 
               ctimeR = timeR - lagDiffR/rec.fs; 
           end
           
           rec.dataL = cdataL;
           rec.dataR = cdataR;

           rec.timeL = ctimeL;
           rec.timeR = ctimeR;
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

            timestamps  = rec.timestamps + rec.offset;
            timestamps2 = timestamps + rec.transition + rec.offset;
            
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
    end % methods
end % recordclass
