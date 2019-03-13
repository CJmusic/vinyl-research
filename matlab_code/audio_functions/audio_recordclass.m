% audio_recordclass
% christopher zaworski
% last edit : march 13 2019 
%
% This is the matlab implementation fo the record class. 
%



classdef audio_recordclass
    properties 
        data = [];
        time = [];
        fs = 0;

        clicks = [];
        path;
    end
    methods 
        function rec = audio_recordclass(file_path);
            path = file_path;
            [data, fs] = audioread(path);
            clicks = audio_clickdetect(data, fs);
        end
        %function r = func(in);
        %    r = 0;
        %end
    end
end
