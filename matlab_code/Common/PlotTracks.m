


function output = PlotTracks(filename);
    [tracks, info_array] = SeperateTracks(filename);
    signal_names = tracks.keys;
    signals = tracks.values;
    mkdir(filename(1:end-4))

    for t = (1:length(tracks))
        trackname = signal_names{t};
        % audiowrite(strcat(filename(1:end-4),'/',signal_names{t},'.wav'),signals{t},96000);  
        y = signals{t};
        x = (0:length(y) - 1)/96000;


        fig = figure('Visible', 'off')
        subplot(2,1,1)
        title(strcat(trackname, ' left channel'))
        plot(x,y(:,1),'k')
        subplot(2,1,2)
        plot(x,y(:,2),'k')
        title(strcat(trackname, ' right channel'))
        grid on;
        ylim([-1,1])
        xlabel('time [s]')
        ylabel('signal level')
        plotname = strcat(filename(1:end-4),'/', signal_names{t}, '.png')
        saveas(fig, plotname);
    
    
    
    
        
    
    
    end


end

