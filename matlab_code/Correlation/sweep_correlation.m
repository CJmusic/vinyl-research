addpath('/Users/cz/Code/vinyl-research/matlab_code/Common/')

sweep_names = {'sweep', 'sweepL', 'sweepR', 'sweepV','sweep2', 'sweepL2', 'sweepR2', 'sweepV2'}

% ref_tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/coherencetest/A0000B0000/031418_A0000B0000r027a1553.770.wav')
% folder = '/Volumes/AUDIOBANK/audio_files/coherencetest/A0000B0000/';

folder = '/Volumes/AUDIOBANK/audio_files/SameRecordTest2';
% ref_tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/SameRecordTest2/r300a1550.844.wav')
ref_tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/SameRecordTest2/r300a1553.621.wav')


files = dir(fullfile(folder,'*.wav'))

keys(ref_tracks)

for i = (1:length(files)) %%loop through records
    filename = files(i).name;
    tracks = SeperateTracks(strcat(files(i).folder,'/',filename));

    for ii = (1:length(sweep_names))
        trackname = sweep_names{ii}

        sweep_ref = ref_tracks(trackname);
        sweep = tracks(trackname);

        if strcmp(trackname, 'sweepR') || strcmp(trackname,'sweepR2');
            sweep_ref = sweep_ref(:,2);
            sweep = sweep(:,2);
        else 
            sweep_ref = sweep_ref(:,1);
            sweep = sweep(:,1);
        end


        [acor_L,lags_L2] = xcorr(sweep_ref,sweep);
        [M_L,I_L] = max(abs(acor_L));
        lagdiff_L2 = lags_L2(I_L);
        lagdiff2 = lagdiff_L2    


        figure(100 + i)
        plot(lags_L2, acor_L)
        [M_L,I_L] = max(abs(acor_L));

        LAGDIFFS(i,ii) = lagdiff2;

    end
end

LAGDIFFS
plot(LAGDIFFS)
