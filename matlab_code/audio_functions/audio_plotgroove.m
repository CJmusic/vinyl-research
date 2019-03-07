
function audio_plotgroove(data)
    clf(figure(1))
    clf(figure(2))
    AUDIO_FILES = {};
    
    %~~~Segmenting Parameters~~~%
    time = time;
    rotation_speed = 33.33333;%45;
    T = 60/rotation_speed; %this is the length of one groove segment
    n_sam = round(T*fs)
    time_seg = time(1:n_sam);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    
    %~~Detect_Signal Parameters~~~%
    winSize= 2^8;%round(fs*dur);
    overlap=0;%round(winSize/2);
    fftsize=2^8;%winSize;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    
    %~~~~~GROOVE PLOTTING~~~~~~%
    seg_array = [];
    name_files= [ {'5.2.wav'},name_files];
    
    figure(2);hold on; legend;
    for i=(1:length(AUDIO_FILES))
        data = AUDIO_FILES{i};
        num_segs = (floor(length(data)/fs/T))
        for ng = 1:num_segs
            %seg_array(:,:,ng) = data(1+(ng-1)*n_sam:ng*n_sam,:);
            data_seg = data(1+(ng-1)*n_sam:ng*n_sam,:);
            [s,f,t] = spectrogram(data_seg,winSize,overlap,fftsize,fs_file,'yaxis');
            lo_freq = sum(s(1:5));
            hi_freq = sum(s(6:end));
            %if abs(hi_freq/lo_freq) < 2.0; %%This is where the DETECT SIGNAL algorithm will go once its active
            if ng > 1 &  ng < 4; %for now just plot grooves 2 & 3  
                plot(time_seg,data_seg,  'DisplayName', [ name_files{i},'groove', num2str(ng)]);
            end
        end    
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~%
end %function plotting
%~~~~~~~~~~~~~~~~~~FUNCTIONS END~~~~~~~~~~~~~~~~~~~~~~%


