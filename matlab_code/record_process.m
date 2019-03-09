% This file will look at a sample recording of our test record and perform all the necessary measurements/analysis. 
% Christopher Zaworski
% Last edit : March 7, 2019
%
%


record_dir = '/Users/cz/OneDrive\ -\ University\ of\ Waterloo/Vinyl_Project/audio_files/MMDDYY_A0000B0000'; 
%this is the directory that the records are recorded under, an example is provided

%process the string and pull info: 
pressing_date = 0; 
recording_date = 0; 
recording_timestamp = 0; 

top_stamper = 0;
top_stamper_hits = 0;
bottom_stamper = 0;
bottom_stamper_hits = 0;

%locate in the csv file

path_csv = 0; 
csv_file = 0; 

%loop through each file: 
for file in record_dir; 
    record_id = [];

    [data, time, fs] = audio_loadfile(file); 
