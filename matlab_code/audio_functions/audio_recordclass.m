% audio_recordclass
% christopher zaworski
% last edit : march 30 2019 
%
% This is the matlab implementation fo the record class. 
%

classdef audio_recordclass < handle %inheriting handle allows methods to update properties
    properties 
        directory = '';
        filename = '';
        extension = '';


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

        timestamps_ref = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
        % this is how many seconds each signal is according to Chris Muth's track listing
        lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
        signal_names = {'leadin','1kHz', '10kHz', '100Hz', 'sweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition', '1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2', '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout'};
        % offset = 0
        offset = 4.25; % as measured on /020818_A0000B0000/02072019_A0000B000r25-A.wav
        % offset = 10.625; % as measured on /020818_A0000B0000/02072019_A0000B000r27-A.wav
        transition = 517.375; % as measured on /020818_A0000B0000/02072019_A0000B000r27-A.wav
        % transition = 518.25; % as measured on /020818_A0000B0000/02072019_A0000B000r25-A.wav
        % lockoutClipped = 953.746;
        lockout = 953.770; % as measured on 031418_A0000B0000r27a %probably the same one

        lagdiff = 0;
        timediff = 0; % this is a time diff calculated as compared to a reference 
        tracks;
        track_times;
    end % properties

    methods 
        function rec = audio_recordclass(file_path);
            disp('record class constructor')
            % strip the filename and path 
            [rec.directory, rec.filename, rec.extension] = fileparts(file_path);
            % strcat(directory,'/',filename,'.mat'); 
            % check if there's already a .mat file of the same name in the directory 
            if exist(strcat(rec.directory,'/',rec.filename,'.mat')) == 10; % 2 is the proper value
                disp('.mat file found, loading....')
                obj = load(strcat(rec.directory,'/',rec.filename,'.mat'));

                rec.data = obj.data;
                rec.dataL = obj.dataL;
                rec.dataR = obj.dataR;
    
                rec.time = obj.time;
                rec.fs = obj.fs;
    
                % rec.clicks = obj.clicks;
                rec.clicksL = obj.clicksL;
                rec.clicksR = obj.clicksR;
    
                rec.click_matrix = obj.click_matrix;
    
                rec.signalsL = obj.signalsL; %
                rec.signalsR = obj.signalsR; %
                rec.signals = obj.signals;
                rec.signal_times = obj.signal_times;
    
                rec.timestamps = obj.timestamps;
                rec.lengths = obj.lengths; 
                rec.signal_names = obj.signal_names; 
    
    
                rec.offset = obj.offset; 
                rec.transition = obj.transition;
                rec.lagdiff = obj.lagdiff;
                rec.timediff = obj.timediff;                 
                rec.tracks = obj.tracks;
                rec.track_times = obj.track_times;
    
            else
            disp('no .mat file found, loading in audio data')

            file_path
            [data, rec.fs] = audioread(file_path);
            rec.data = data;
            rec.time = (0:length(rec.data)-1)/rec.fs;
            % rec.dataL = data(:,1);
            % rec.dataR = data(:,2);
            %rec.clicksL = audio_clickdetect(rec.dataL, rec.fs);
            %rec.clicksR = audio_clickdetect(rec.dataR, rec.fs);

            end % if .mat file exists
        end % constructor

        function obj = save_obj(rec);
            disp('inside save_obj')
            obj.data = rec.data;
            obj.dataL = rec.dataL;
            obj.dataR = rec.dataR;

            obj.time = rec.time;
            obj.fs = rec.fs;

            obj.clicks = rec.clicks;
            obj.clicksL = rec.clicksL;
            obj.clicksR = rec.clicksR;

            obj.click_matrix = rec.click_matrix;

            obj.signalsL = rec.signalsL; %
            obj.signalsR = rec.signalsR; %
            obj.signals = rec.signals;
            obj.signal_times = rec.signal_times;

            obj.timestamps = rec.timestamps;
            obj.lengths = rec.lengths; 
            obj.signal_names = rec.signal_names; 


            obj.offset = rec.offset; 
            obj.transition = rec.transition;
            obj.lagdiff = rec.lagdiff;
            obj.timediff = rec.timediff; 
            obj.tracks = rec.tracks;
            obj.track_times = rec.track_times;

            % type(obj)
            obj            
            path = strcat(rec.directory,'/', rec.filename,'.mat')
            save(path, '-v7.3', '-struct','obj')
        end % save_recordclass

        function lagcorrect(rec, ref);
        % function lagcorrect(rec); %%% WORKING %%%%
            % this function should correct the record for any lagdiff calculated by another method: 
            % the method currently implemented is okay, it simply adds the lag 
            % to the track times, need to reslice the tracks based on it

            %%% WORKING $$$
            % for k = (keys(rec.track_times));
            %     % rec.track_times(i)
            %     rec.track_times(k{1}) = rec.track_times(k{1}) + rec.lagdiff/rec.fs;
            % end
            %%% WORKING ENDS %%% 
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
            rec.track_times = {};
            rec.signals = {};
            rec.signal_times = {};
            timestamps = [];
            timestamps2 = [];
            signalsL = [];
            signalsR = [];
            time_seg = [];


            signal_names = {'leadin','1kHz', '10kHz', '100Hz', 'sweep', 'quiet', '3150Hz', '1kHzL', 'sweepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV','transition', '1kHz2', '10kHz2', '100Hz2', 'freqsweep2', 'quiet2', '3150Hz2', '1kHzL2', 'sweepL2', '1kHzR2', 'sweepR2', '1kHzV2', 'sweepV2','leadout'};
            % timestamps  = rec.timestamps + rec.offset + rec.timediff;
            % timestamps2 = timestamps + rec.transition + rec.offset + rec.timediff;
            % disp('timestamps')
            disp('timestamps')
            timestamps  = [rec.timestamps_ref + rec.offset + rec.lagdiff/rec.fs, rec.timestamps_ref + rec.offset + rec.transition + rec.lagdiff/rec.fs] %+ rec.timediff;
            timestamps
            signal_times(1) = [1,floor(timestamps(1)*rec.fs)]

            for i=(2:length(timestamps)-1)
                signal_times(i) = [floor(timestamps(i-1))*rec.fs,floor(timestamps(i))*rec.fs];
            end

            signal_times(length(timestamps)) = [floor(timestamps(length(timestamps)))*rec.fs,length(rec.data)];
            signal_times
            % % get the needledrop
            % % signalsL = rec.dataL(1:floor((timestamps(1))*rec.fs));
            % % signalsR = rec.dataR(1:floor((timestamps(1))*rec.fs));

            % time_seg = rec.time(1:floor((timestamps(1))*rec.fs));
            
            % % rec.signals{end + 1} = [signalsL, signalsR];
            
            % rec.signal_times{end + 1} = [time_seg];
            
            % % first set of signals on the disk
            % for i = (1:length(timestamps)-1);            

            % end
            
            % % the extended silence section 'transition' between the two sets of signals 
            % % disp('transition track')
            % % rec.transition
            % % timestamps(end)
            % % timestamps2(1)
            
            % % signalsL = rec.dataL(floor(timestamps(end)*rec.fs):floor(timestamps2(1)*rec.fs));
            % % signalsR = rec.dataR(floor(timestamps(end)*rec.fs):floor(timestamps2(1)*rec.fs));
            
            % time_seg = rec.time(floor(timestamps(end)*rec.fs):floor((timestamps2(1))*rec.fs));
	
            % rec.signals{end + 1} = [signalsL, signalsR];
            % rec.signal_times{end + 1} = [time_seg];
            
            % % % second set of signals on the disk
            % % for i = (1:length(timestamps2)-1);
            % %     signalsL = rec.dataL(floor(timestamps2(i)*rec.fs):floor(timestamps2(i+1)*rec.fs));
            % %     signalsR = rec.dataR(floor(timestamps2(i)*rec.fs):floor(timestamps2(i+1)*rec.fs));
            % %     time_seg = rec.time(floor(timestamps2(i)*rec.fs):floor(timestamps2(i+1)*rec.fs));
            % %     rec.signals{end + 1} = [signalsL, signalsR];
            % %     rec.signal_times{end + 1} = [time_seg];
            % % end
            
            % % % get the rest of the file
            % % signalsL =rec.dataL(floor(timestamps2(end)*rec.fs):end);
            % % signalsR = rec.dataR(floor(timestamps2(end)*rec.fs):end);
            % % time_seg = rec.time(floor(timestamps2(end)*rec.fs):end);
            % % rec.signals{end + 1} = [signalsL, signalsR];
            % rec.signal_times{end + 1} = [time_seg];
           
            % % disp('Size of signal names, then signals, then times')
            % % size(rec.signal_names)
            % % size(rec.signals)

            % rec.tracks = containers.Map(rec.signal_names, rec.signals);
            % rec.track_times = containers.Map(rec.signal_names, rec.signal_times);
            % % disp('printing tracks')
            
            disp('done processing tracks')

        end % function signals

        function leadout_track(rec);

        end


        function transition_track(rec);
            % this function grabs the approximate position of the transition track for the purpose of 
            % lining up the file to a reference
            % rec.tracks = 0; % clear any previous tracks info
            % rec.signals = {};
            
            timestamps  = rec.timestamps_ref + rec.offset + rec.timediff;
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
            % rec.clicks = audio_clickdetect(rec.data, rec.fs); 
            rec.clicksL = audio_clickdetect(rec.data(:,1), rec.fs); 
            rec.clicksR = audio_clickdetect(rec.data(:,2), rec.fs); 
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
