function[rms_value]=rms_response(array)
% [rms_value]=rms_response(array)
% -------------------------------
% This function computes the root-mean-square value of an array.
% See also the function 'avg_response'.
% John Vanderkooy June 2019
%-------------------------------------------------
numel=length(array);
rms_value=sqrt(sum(abs(array).^2)/numel);


