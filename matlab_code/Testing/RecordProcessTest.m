close all; clear all; clc;


disp('~~~~~~~~~~~~TESTING RECORDPROCESS~~~~~~~~~~~~')

addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\')
addpath('D:\Code\vinyl-research\matlab_code\audio_functions')
file = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\03141_A0000B0000r030b.wav'

recordProcess(file)
