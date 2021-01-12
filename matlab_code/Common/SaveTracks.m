


function output = SaveTracks(filename);
    [tracks, info_array] = SeperateTracks(filename);
    signal_names = tracks.keys;
    signals = tracks.values;
    mkdir(filename(1:end-4))

    for t = (1:length(tracks))
        audiowrite(strcat(filename(1:end-4),'/',signal_names{t},'.wav'),signals{t},96000);
    end

end