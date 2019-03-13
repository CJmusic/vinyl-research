% audio_recordclass
% christopher zaworski
% last edit : march 13 2019 
%
% This is the matlab implementation fo the record class. 
%



classdef audio_recordclass
    properties 
        dataL = [];
        dataR = [];
        time = [];
        fs = 0;

        clicksL = [];
        clicksR = [];

        signalsL = []; %
        signalsR = []; %
        
            path;
    end % properties
    methods 
        function rec = audio_recordclass(file_path);
            path = file_path;
            [data, rec.fs] = audioread(path);
            rec.dataL = data(:,1);
            rec.dataR = data(:,2);
            rec.clicksL = audio_clickdetect(rec.dataL, rec.fs);
            rec.clicksR = audio_clickdetect(rec.dataR, rec.fs);
        end
        function signals  = signal_process();
            signals = {'leadin','1kHz', '10kHz', '100Hz', 'freqsweep', 'quiet', '3150Hz', '1kHzL', 'swepL', '1kHzR', 'sweepR', '1kHzV', 'sweepV', 'leadout'} 
            timesstamps = []; %need detect vinyl noise, probably with detect signal  
             
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
        
        end
    end % methods
end % recordclass
