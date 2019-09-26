function [avg_value]=avg_response(array)
% [avg_value]=avg_response(array)
% -------------------------------
% This function computes the average of the magnitude values of an array.
% A sinusoidal signal of amplitude A will give a value A/sqrt(2), like rms.
% See also the function 'rms_response'.
% John Vanderkooy June 2019
%-------------------------------------------------
numel=length(array);
avg_value=(pi/(2*sqrt(2)))*sum(abs(array))/numel;


