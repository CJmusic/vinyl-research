
clear all; clc;close all;
disp('----------------start of program--------------------')
set(0,'DefaultLineLinewidth',1.5)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
%
try
    pkg load signal %for Octave
catch
end

addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')