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

        time = [];
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
        signal_times = {};
        
        %%% RECORD INFO 
        % timestamps = [0, 2, 62, 92, 124, 160, 182, 248, 268, 306, 326, 364, 384, 419.5];
        timestamps = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
        % this is how many seconds each signal is according to Chris Muth's track listing
        lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
        signal_names = {'leadin','1kHz', '10kHz', '100Hz', 'freqsweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition', '1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2', '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout'};

        % offset = 4.25; % as measured on /020818_A0000B0000/02072019_A0000B000r25-A.wav
        offset = 10.625; % as measured on /020818_A0000B0000/02072019_A0000B000r27-A.wav
        transition = 517.375; % as measured on /020818_A0000B0000/02072019_A0000B000r27-A.wav
        % transition = 518.25; % as measured on /020818_A0000B0000/02072019_A0000B000r25-A.wav
        lagdiff = 0;
        timediff = 0; % this is a time diff calculated as compared to a reference 
        tracks;
        track_times;
    end % properties
    methods 
        function rec = audio_recordclass(file_path);
            disp('record class constructor')
            file_path
            [data, rec.fs] = audioread(file_path);
            rec.data = data;
            rec.time = (0:length(rec.data)-1)/rec.fs;
            rec.dataL = data(:,1);
            rec.dataR = data(:,2);
            %rec.clicksL = audio_clickdetect(rec.dataL, rec.fs);
            %rec.clicksR = audio_clickdetect(rec.dataR, rec.fs);
        end

        function lagcorrect(rec);
            % this function should correct the record for any lagDiff calculated by another method: 
            % the method currently implemented is okay, I want to try circshift
            % also see the matlab function lag !!!! 
            disp('inside lagcorrect')
            rec.lagdiff

            if rec.lagdiff == 0; %then nothing needs to be corrected  
                disp('no lag diff')
            else;
                rec.data = circshift(rec.data, -1*rec.lagdiff);
                rec.dataL = rec.data(:,1);
                rec.dataR = rec.data(:,2);
            end
            
            rec.time = (0:length(rec.data)-1)/rec.fs;
            % time = (0:length(rec.data)); 
        %    if lagDiff > 0;
        %        % positive lagDiffL means that the data array is delayed compared to the reference
        %        cdata = rec.data(lagDiff + 1:end,:); 
        %        ctime = time - lagDiff/rec.fs; 
        %    elseif lagDiff < 0; 
        %        % negative lagDiffL means that the data array is ahead of the reference  
        %        cdataL = rec.data(1: end + lagDiff,:);   
        %        ctimeL = time - lagDiff/rec.fs;
        %    elseif lagDiff == 0; %then nothing needs to be corrected  
        %        cdata = rec.data; %no lag 
        %        ctime = time - lagDiff/rec.fs; 
        %    end

           
        %    rec.data = cdata;

        %    rec.time = ctime;
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

            signal_names = {'leadin','1kHz', '10kHz', '100Hz', 'freqsweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition', '1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2', '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout'};
            % timestamps  = rec.timestamps + rec.offset + rec.timediff;
            % timestamps2 = timestamps + rec.transition + rec.offset + rec.timediff;
            disp('timestamps')
            timestamps  = rec.timestamps + rec.offset + rec.timediff
            timestamps2 = rec.timestamps + rec.offset + rec.transition + rec.timediff
            
            % get the needledrop
            signalsL = rec.dataL(1:floor((timestamps(1))*rec.fs));
            signalsR = rec.dataR(1:floor((timestamps(1))*rec.fs));
            time_seg = rec.time(1:floor((timestamps(1))*rec.fs));
            rec.signals{end + 1} = [signalsL, signalsR];
            rec.signal_times{end + 1} = [time_seg];
            
            % first set of signals on the disk
            for i = (1:length(timestamps)-1);
                disp('FOR LOOP LENGTHS')
                starting = floor(timestamps(i))
                ending = floor(timestamps(i+1))
                lengthing =  ending - starting 
                signalsL = rec.dataL(floor(timestamps(i)*rec.fs):floor(timestamps(i+1)*rec.fs));
                signalsR = rec.dataR(floor(timestamps(i)*rec.fs):floor(timestamps(i+1)*rec.fs)); 
                time_seg = rec.time(floor(timestamps(i)*rec.fs):floor(timestamps(i+1)*rec.fs));
                rec.signals{end + 1} = [signalsL, signalsR];
                rec.signal_times{end + 1} = [time_seg];
            end
            
            % the extended silence section 'transition' between the two sets of signals 
            disp('transition track')
            rec.transition
            timestamps(end)
            timestamps2(1)
            signalsL = rec.dataL(floor(timestamps(end)*rec.fs):floor(timestamps2(1)*rec.fs));
            signalsR = rec.dataR(floor(timestamps(end)*rec.fs):floor(timestamps2(1)*rec.fs));
            time_seg = rec.time(floor(timestamps(end)*rec.fs):floor((timestamps2(1))*rec.fs));
	
            rec.signals{end + 1} = [signalsL, signalsR];
            rec.signal_times{end + 1} = [time_seg];
            
            % second set of signals on the disk
            for i = (1:length(timestamps2)-1);
                signalsL = rec.dataL(floor(timestamps2(i)*rec.fs):floor(timestamps2(i+1)*rec.fs));
                signalsR = rec.dataR(floor(timestamps2(i)*rec.fs):floor(timestamps2(i+1)*rec.fs));
                time_seg = rec.time(floor(timestamps2(i)*rec.fs):floor(timestamps2(i+1)*rec.fs));
                rec.signals{end + 1} = [signalsL, signalsR];
                rec.signal_times{end + 1} = [time_seg];
            end
            
            % get the rest of the file
            signalsL =rec.dataL(floor(timestamps2(end)*rec.fs):end);
            signalsR = rec.dataR(floor(timestamps2(end)*rec.fs):end);
            time_seg = rec.time(floor(timestamps2(end)*rec.fs):end);
            rec.signals{end + 1} = [signalsL, signalsR];
            rec.signal_times{end + 1} = [time_seg];
           
            disp('Size of signal names, then signals')
            size(rec.signal_names)
            size(rec.signals)
            rec.tracks = containers.Map(rec.signal_names, rec.signals);
            rec.track_times = containers.Map(rec.signal_names, rec.signal_times);
            disp('printing tracks')
            disp('done processing tracks')

        end % function signals
        function transition_track(rec);
            % this function grabs the approximate position of the transition track for the purpose of 
            % lining up the file to a reference
            % rec.tracks = 0; % clear any previous tracks info
            % rec.signals = {};
            
            timestamps  = rec.timestamps + rec.offset + rec.timediff;
            timestamps2 = timestamps + rec.transition + rec.offset + rec.timediff;


            % the extended silence section 'transition' between the two sets of signals 
            signalsL = rec.dataL(floor(timestamps(end)*rec.fs): ...
                                           floor((timestamps2(1))*rec.fs));
            
            signalsR = rec.dataR(floor(timestamps(end)*rec.fs): ...
                                           floor((timestamps2(1))*rec.fs));
            transition = [signalsL, signalsR];
            disp('inside rec class')
            size(transition)
            rec.tracks = containers.Map('transition', transition);
            class(rec.tracks)
            size(rec.tracks('transition'))
        end % transition track
        
        function clickdetect(rec);
            rec.clicks = audio_clickdetect(rec.data, rec.fs); 
        end % function detect_clicks

        function clicklineup(rec, clicks_ref)
            [rec.click_matrix, rec.lagdiff] = audio_clickmatrix(rec.clicks, clicks_ref);
            disp('lagdiff')
            size(rec.click_matrix)
            % size(rec.lagdiff)
            rec.lagdiff
            rec.timediff = rec.lagdiff/rec.fs;
            % rec.data = circshift(rec.data, rec.lagdiff);
            
        end % function click_lineup
    end % methods
end % recordclass
