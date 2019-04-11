% audio_rec.ecordclass.m 
% christopher zaworski
% last edit: march 17 2018 
%
% This is the implementation of the record class. It is a
% special case of a record class used to line up other audio files. 
%

% This class is currently tuned to use the following: 
% /020818_A0000B0000/02072019_A0000B000r27-A.wav
% as the file

classdef audio_refrecordclass < audio_recordclass
    properties 
           % signals = {};

           % %%% RECORD INFO 
           % timestamps = [0, 2, 62, 92, 124, 160, 182, 248, 268, 306, 326, 364, 384, 419.5];
           % % this is how many seconds each signal is according to Chris Muth's track listing
           % lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
           % signal_names = {'needledrop','leadin','1kHz', '10kHz', '100Hz', 'freqsweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition','leadin2','1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2', '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout2'};

           % offset = 4.25; % as measured on /020818_A0000B0000/02072019_A0000B000r25-A.wav
           % %time = (0:length(dataL)-1)/rec.fs;
           % transition = 518.25; % as measured on /020818_A0000B0000/02072019_A0000B000r25-A.wav
           % tracks = 0;
    end % properties

    methods
        function rec = audio_refrecordclass(file_path); 
            disp('ref recordclass constructor')
            file_path
            %call the record class constructor 
            rec@audio_recordclass(file_path);
        %end
        %function tracks = track_listings(); 
            %rec = caudio_recordclass(rec.file_path);

            timestamps  = rec.timestamps + rec.offset;
            timestamps2 = timestamps + rec.transition;
            
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
        end % constructor
    end % methods
end % classdef
