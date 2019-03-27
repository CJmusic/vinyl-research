%{

This file can also count the clicks in the record and see what clicks are common 
across all records (ones that are in the reference file as well), and ones that 
aren't. 

%}
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/260219_noisereferenceinst/');

addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/');
audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_A0000B0000r26fivetrials/'

AUDIO_FILES = {'one.wav','two.wav','three.wav', 'four.wav', 'five.wav'};
%AUDIO_FILES = {[fsinst_shorted, fsinst_shortedgained, gain10reference, gain10recordnoise, recordnoise, reference, system_noise, system_gained] };

path_ref = strcat(audio_dir,AUDIO_FILES{2}); % choose a file to use as the reference
path_ref
[data_ref, fs_ref] = audioread(path_ref);
time_ref = (0:length(data_ref)-1)/fs_ref;

t_s = 4%67.0;% seconds, the start of the region in the file being analyzed
t_e = 30%70.0;% seconds, the end of the region being analyzed

data_ref = data_ref(t_s*fs_ref:t_e*fs_ref,1);
time_ref = (0:length(data_ref)-1)/fs_ref; 

[ clicks_ref ] = audio_clickdetect(data_ref, fs_ref);

clf(figure(1)); clf(figure(2)); %clear any figures used for plotting 

for i = (1:2);%length(AUDIO_FILES));         
    [data, time, fs] = audio_load(strcat(audio_dir,AUDIO_FILES{i}));
    data = data(t_s*fs:t_e*fs,1);
    time = (1:length(data))/fs;
    %[data, time] = audio_lineup(data, fs, time, data_ref);
    [ clicks ] = audio_clickdetect(data, fs);
    %lagdiff = mode(clicks)

    time = (1:length(data))/fs;
       %[clicks_ref] = audio_clickdetect(data_ref, fs_ref);

    diff_array = [] % this array contains the distances between every click, each row 
                     % represents a click in the referenc each column represents a click in the file being looked at 
    for xi = (1:length(clicks_ref));%this makes an array with the distance between each click in the two
                                    % files
        diff_array = [diff_array; clicks - clicks_ref(xi)];
    end

    lagdiff = mode(diff_array(:)) % the time difference between the two signals is the most common distance between clicks
    timediff = lagdiff/fs


    if lagdiff > 0;
        % positive lagdiff means that the data array is delayed compared to the reference
        cdata = data(lagdiff + 1:end,:); 
        time = time - lagdiff/fs; 
    elseif lagdiff < 0; 
        % negative lagdiff means that the data array is ahead of the reference  
        cdata = data(1: end + lagdiff,:);   
        time = time - lagdiff/fs;
    elseif lagdiff == 0; %then nothing needs to be corrected  
        cdata = data; %no lag 
        time = time - lagdiff/fs; 
    end

    size_diff = length(data_ref) - length(cdata)
    cdata_ref = data_ref(abs(size_diff)+1:end,:); %need to      
    %if size_diff > 0; %% the following is to ensure the data arrays are the same length for calculating the coherence
    %                  %% probably this isn't necessary after the lag diff is properly handled 
    %    disp('>0'); 
    %    cdata_ref = data_ref(size_diff+1:end,:); 
    %elseif size_diff < 0; 
    %    cdata = cdata(size_diff,:);
    %    cdata_ref = data_reff;
    %else;
    %    cdata_ref = data_ref;
    %end

    disp('sizes corrected')
    size(cdata)
    size(cdata_ref)

    [amp_coh, freq_coh] = audio_mscohere(cdata_ref, cdata, fs);

    % plotting below : 
    figure(1);
    grid on; hold on;
    plot(time,data(:,1));%, 'Color', [i/5,i/5,i/5] );
    legend(AUDIO_FILES);
    title('Click detection: Amplitude vs. Time'); 
    xlabel('Time [s]');
    ylabel('Amplitude');

    figure(2); grid on; hold on;
    %plot(freq_coh, amp_coh); 
    %set(gca, 'XScale', 'log');
    semilogx(freq_coh, amp_coh);
    legend(AUDIO_FILES);
    xlabel('frequency [Hz]')
    title('Coherence compared to first recording')
end


