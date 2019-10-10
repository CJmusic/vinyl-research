
addpath('E:\audio_files\A0000B0000\')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/') 
% audioFile = ('E:\audio_files\A0000B0000\A0000B0000.csv')
% pressFile = ('E:\audio_files\A0000B0000\oct10A0000B0000.csv')
folder = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/')
rawFile = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/A0000B0000.csv')
timeFile = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/oct10A0000B0000.csv')


% read in the press file and the time file 
% join the two arrays 


%#####_SensorValues
% tracks the data from the sensors in the press via timestamps


%#####_JobDetailsCurrent
% tracks cycle times


%#####_ADAPT_DATA
% tracks the press options via timestamps

%#####_TimeStamps
% Recorded by hand timestamps


% get the hand recorded timestamp
% identify a pressing cycle by press position (one open and close)
% -> need to keep in mind times when parameter's changed 



rawData = readtable(pressFile);
timeData = readtable(timeFile);

