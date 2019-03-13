% audio_detectsignal.m
% christopher zaworski
% last edit : march 13 2018
%
% this file takes a data array and outputs an array with 0 
% where there is record silence, 1 where there are signals


function signal_array = audio_detectsignal(data, fs); 
    signal_array = [];

    % wait until a sample reaches above that value, then wait a few samples to 
    % see of it another sample comes along that reaches that value, as that happens put 
    % signal there

    threshold = 0.2; % simple threshold for now

    for i = (1:length(data));
        %i
        %length(signal_array)
        if data(i) > threshold;
            signal_array = [signal_array, 1];
        %elseif i > 65 && sum(signal_array(i-65:i-1)) > 1 ;
            %signal_array = [signal_array, 1];
        else 
            signal_array = [signal_array, 0]; 
        end
    end
    clf(figure(1));
    figure(1);hold on; grid on;
    plot(data);
    plot(signal_array,'r');
end

