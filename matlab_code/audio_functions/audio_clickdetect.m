%{
Lining up clicks: 

clicks = [ sample numbers ] 

% a click has a sustained high derivative 
% the energy of a click is between 2-4 kHz

%} 

%% TESTING THE CLICK DETECT FUNCTION
%addpath('audio_functions')
%addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/260219_noisereferenceinst/');
%addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/');
%audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/'
%
%AUDIO_FILES = {'one.wav','two.wav','three.wav', 'four.wav', 'five.wav'};

% clear all; clc; close all;
% addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/');
% audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/';

% audio_bin = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r26-96kHz.wav';

% AUDIO_FILES = {'Bcorrelation_test_1.wav','Bcorrelation_test_2.wav','Bcorrelation_test_3.wav'};

% %[data, time, fs] = audio_load(strcat(audio_dir,AUDIO_FILES{1}));
% [data, time, fs] = audio_load(audio_bin);
% data = data(5.1*fs : 8.1*fs,:);
% audio_clickdetecttest(data, fs);

% function clicks = audio_clickdetecttest(data, fs);
%     % implement a band pass filter 1-2 kHz
%     %fc = 1000.0 
%     %  [A,B,C,D] = butter(10,[1000 10000]/fs/2);
%     %  d = designfilt('bandpassiir','FilterOrder',20, ...
%     %    'HalfPowerFrequency1',1000,'HalfPowerFrequency2',10000, ...
%     %    'SampleRate',fs);
%     %
%     %  sos = ss2sos(A,B,C,D);
%     %  fvt = fvtool(sos,d,'Fs',1500);
%     %  legend(fvt,'butter','designfilt')
%     time = (1:length(data))/fs;
%     [b, a] = butter(1, [1000.0 4000.0]/(fs/2),'bandpass'); % this should 
%     data_f = filter(b, a, data); % use filtfilt
%     figure(1);
%     freqs(b,a)
%     title('bandpass filter')
%     figure(2); grid on; hold on
%     plot(time, data_f)
%     title('audio data')
%     % find the start of the click 
%     p_data = data_f.^2;
%     figure(3); grid on; hold on;
%     plot(time, p_data)
%     title('audio power')
 
%     clicks = [];
%     zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);% Returns Zero-Crossing Indices Of Argument Vector
%     zero_indexes = zci(data);
%     size(zero_indexes)
%     size(p_data)
%     avg_peak = 0.0;
%     prev_peak = 1.0;
%     peak_value = 0.0;
%     %[peak_values, peak_indices] = max(p_data(zero_indexes));
%     peak_values = [];   
%     peak_indices = [];
%     disp('peak values')
%     size(peak_values)
%     %[peak_values, clicks] = findpeaks(p_data(:,1));
%     for i = (1:length(zero_indexes)-1);
%         [peak_value, peak_index] = max(p_data(zero_indexes(i): zero_indexes(i+1)));
%         peak_values = [peak_values, peak_value];
%         peak_indices = [peak_indices, peak_index];
%         %if peak_value > 10.0*prev_peak;
%        % if peak_value/prev_peak > 10.0;
%        %    click = zero_indexes(i) + peak_index;
%        %    clicks = [clicks, click];
%        % else; 
%        %    avg_peak = (avg_peak + prev_peak)/2;
%        % end
%        % prev_peak = peak_value;
%     end
%     size(peak_values)
%     %avg_peak = mean(peak_value) 
%     %clicks
%     size(clicks)
%     %  plot(time(1:end-t_diff),d_data); grid on; hold on;
%     figure(2); hold on;
%  %   for xi = 1:length(clicks)-1;
%  %       x1 = time(clicks(xi));
%  %       line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
%  %   end
%  %   figure(3); hold on;
%  %   for xi = 1:length(clicks)-1;
%  %       x1 = time(clicks(xi));
%  %       line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
%  %   end
%     % determine the click's polarity 
    
%     % length and number of Oscillations 
    
% end   


% CODE BELOW is old click method of first differences  


function clicks = audio_clickdetect(data, fs);
    time = (0:length(data)-1)/fs; 
    d_data = diff(data,1)*fs;% since delta_t = 1/fs:w;
    clicks = [];
    threshold = 10*rms(d_data(1:1025));
    disp('starting click detect')
    for i = (1:length(data)-1);
        if i > 512 && i + 512 < length(d_data);
            threshold = threshold - abs(d_data(i-511)) + abs(d_data(i+511));
            % threshold = 10*rms(d_data(i-511:i+512))
        end
        if d_data(i) > threshold;
            click = i;
            %click_timestamp = i/fs;   
            clicks = [clicks, click];   
        end
    end
    disp('Number of clicks');
    disp('new file')
    size(clicks) 
    
    %plotting below
    
    %clf(figure(1));clf(figure(2)); 
    %figure(1); grid on; hold on;
    %plot(time, data); 
    % for xi = 1:length(clicks);
    %      x1 = time(clicks(xi));
    %      line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    % end
    % title('Click detection: Amplitude vs. Time'); 
    % xlabel('Time [s]');
    % ylabel('Amplitude');
    
    %figure(2); grid on; hold on;
    %t_diff = length(time) - length(d_data)
    %plot(time(1:end-t_diff),d_data); grid on; hold on;
    %for xi = 1:length(clicks);
    %     x1 = time(clicks(xi));
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    % end
    
    % title('Click detection: 1st Derivative vs. Time'); 
    % xlabel('Time [s]');
    % ylabel('dA/dt (s^-1)');
end 




%}