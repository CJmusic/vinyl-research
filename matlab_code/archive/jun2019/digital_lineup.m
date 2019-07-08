
clc;close all;
disp('-----------digital_lineup.m------------')
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');

path_digital_files = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/reference_files/digital_signals-resampled/'; 

sweep_file = 'sweep16kHz.wav';
tone100_file = 'tone100.wav';
tone1000_file = 'tone1000.wav'; 
tone3150_file = 'tone3150.wav';
tone10000_file = 'tone10000.wav';

%%{
%
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
%
%%}

% name_files = {tone1000_file, tone10000_file, tone100_file, sweep_file  tone3150_file, tone1000_file  sweep_file, tone1000_file, sweep_file, tone1000_file, sweep_file};


%now reference is recording
reference = audio_recordclass('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r28a.wav');
reference.process_tracks();

tone1000 = audio_recordclass(strcat(path_digital_files,tone1000_file));

lagdiff = audio_lineup(reference.tracks('1kHz'),tone1000.data, reference.fs)

time_tone = (0:length(tone1000.data)-1)/tone1000.fs; %- lagdiff/tone1000.fs;

figure(1); grid on; hold on;
plot(reference.track_times('1kHz'),reference.tracks('1kHz'))
plot(time_tone, tone1000.data, 'g');



