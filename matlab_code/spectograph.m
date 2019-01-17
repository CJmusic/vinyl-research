% pkg load signal

[data, fs] = audioread('../../audio_files/300818_TestPressing_t9/b-sides/110918UW_regrindmouldsheated-r6B.wav');
data = data(1:fs*10);
% Window Size = 256, Window Overlap = Default, Frequency Range = Default;
% spectrogram(data, windowSize, windowOverlap, freqRange, fs, 'yaxis');

%Larger Window Size value increases frequency resolution
%Smaller Window Size value increases time resolution
%Specify a Frequency Range to be calculated for using the Goertzel function
%Specify which axis to put frequency

%%EXAMPLE
figure();
[S, f, t] = specgram(data, 256, fs);
imagesc (t, f, log(S));
%Feel free to experiment with these values
