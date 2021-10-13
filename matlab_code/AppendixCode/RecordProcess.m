

function output = recordProcess(file)
            output = [];
            fs = 96000;
            [tracks, info_array] = SeperateTracks(file);
            
            info_array
            lagdiff = info_array(1);
            normalization_L = info_array(2);
            normalization_R = info_array(3);

            signal_names = tracks.keys;
            signals = tracks.values;

            display('LOADING REFERENCE')
            if ismac() == true
                if file(length(file)-4) == 'a'
                    reference = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
                    %% Reference 02072019_A0000B000r27a.wav 
                    offset = 15; 
                elseif file(length(file)-4) == 'b'
                    disp('PC')
                    reference = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028b.wav'); 
                    %% Reference 02072019_A0000B000r27b.wav 
                    offset = 13.1;
                else 
                    disp('NO SIDE FOUND, USING SIDE A REFERENCE')
                    reference = ('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
                    %% Reference 02072019_A0000B000r27a.wav 
                    offset = 15; 
                end
            end
            if ispc() == true
                disp('IS PC')
                if file(length(file)-4) == 'a'
                    reference = ('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
                    %% Reference 02072019_A0000B000r27a.wav 
                    offset = 15; 
                elseif file(length(file)-4) == 'b'
                    disp('PC')
                    reference = ('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028b.wav'); 
                    %% Reference 02072019_A0000B000r27b.wav 
                    offset = 13.1;
                else 
                    disp('NO SIDE FOUND, USING SIDE A REFERENCE')
                    reference = ('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 
                    %% Reference 02072019_A0000B000r27a.wav 
                    offset = 15; 
                end
            end

            
            for t = 1:length(signal_names)
                track_name = signal_names{t};
                sig = signals{t};
                disp(strcat('track  ...',track_name))
    
    
                csig = [];
                CLICKS_R = [];
                CLICKS_L = [];
                RMS_L = [];
                RMS_R = [];
                THD_L = [];
                THD_R = [];
                
                % if t == 1
                [csig(:,1), CLICKS_L] = ClickDetect(sig(:,1));
                [csig(:,2), CLICKS_R] = ClickDetect(sig(:,2));

                REFSa1_L = load('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/clicks/A0000B0000r028a1558.066clicks_L.mat');

                REFSa1_R = load('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/clicks/A0000B0000r028a1558.066clicks_R.mat');

                refsa1_L = REFSa1_L.clicks_L(signal_names{t});
                refsa1_R = REFSa1_R.clicks_R(signal_names{t});

                REFSb1_L = load('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/clicks/A0000B0000r028a1558.066clicks_L.mat');
                REFSb1_R = load('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/clicks/A0000B0000r028a1558.066clicks_R.mat');

                refsb1_L = REFSb1_L.clicks_L(signal_names{t});
                refsb1_R = REFSb1_R.clicks_R(signal_names{t});


                REFSa2_L = load('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/clicks/131a1552.480clicks_L.mat');
                REFSa2_R = load('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/clicks/131a1552.480clicks_R.mat');
                
                refsa2_L = REFSa2_L.clicks_L(signal_names{t});
                refsa2_R = REFSa2_R.clicks_R(signal_names{t});


                REFSb2_L = load('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/clicks/075b1559.017clicks_L.mat');
                REFSb2_R = load('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/clicks/075b1559.017clicks_R.mat');

                refsb2_L = REFSb2_L.clicks_L(signal_names{t});


                refsb2_R = REFSb2_R.clicks_R(signal_names{t});
                
    
                % need to do the reference here by track 
                commonclicksa1_L = CommonClicks(CLICKS_L, refsa1_L);
                commonclicksa1_R = CommonClicks(CLICKS_R, refsa1_R);
                commonclicksb1_L = CommonClicks(CLICKS_L, refsb1_L);
                commonclicksb1_R = CommonClicks(CLICKS_R, refsb1_R);

                commonclicksa2_L = CommonClicks(CLICKS_L, refsa2_L);
                commonclicksa2_R = CommonClicks(CLICKS_R, refsa2_R);
                commonclicksb2_L = CommonClicks(CLICKS_L, refsb2_L);
                commonclicksb2_R = CommonClicks(CLICKS_R, refsb2_R);

    
                
                clicks_L = length(CLICKS_L);
                clicks_R = length(CLICKS_R); 

                rmssig_L = rms(sig(:,1));
                rmssig_R = rms(sig(:,2));
                rmscsig_L = rms(csig(:,1));
                rmscsig_R = rms(csig(:,2));

                RMSclicks_L = 20.0*log10(sqrt(rmssig_L^2-rmscsig_L^2));
                RMSclicks_R = 20.0*log10(sqrt(rmssig_R^2-rmscsig_R^2));

    
                RMS_L = 20.0*log10(rms(csig(:,1)));
                RMS_R = 20.0*log10(rms(csig(:,2)));
    
                Aw = audio_Aweighting(csig);
                CCIRw = audio_CCIRweighting(csig);
    
                A_L = 20.0*log10(rms_response(Aw(:,1)));
                A_R = 20.0*log10(rms_response(Aw(:,2)));
    
                CCIR_L = 20.0*log10(avg_response(CCIRw(:,1)));
                CCIR_R = 20.0*log10(avg_response(CCIRw(:,2)));
    
    
                THD_L = thd(csig(:,1),fs);
                THD_R = thd(csig(:,2),fs);
    
    
                if ismember(signal_names(t), {'1kHzL', 'sweepL', '1kHzL2', 'sweepL2'})             
                    stereo_bleed = StereoBleed(csig,1);
                elseif ismember(signal_names(t), {'1kHzR', 'sweepR', '1kHzR2', 'sweepR2'}) 
                    stereo_bleed = StereoBleed(csig,2);
                else
                    stereo_bleed = 0;
                end
    
                if ismember(signal_names(t), {'3150Hz', '3150Hz2', '1kHz', '1kHz2', '1kHzL', '1kHzL2','1kHzR', '1kHzR2', '1kHzV', '1kHzV2'})
                
                    
                    [test_freq_L, wfreqspecamplitude_L, freqrms_L, WFrms_L] = WowFlutter(csig(:,1));
                    [test_freq_R, wfreqspecamplitude_R, freqrms_R, WFrms_R] = WowFlutter(csig(:,2));

             
                else 
                    test_freq_L = 0;
                    test_freq_R = 0;
                    wfreqspecamplitude_L = 0;
                    wfreqspecamplitude_R = 0;
                    freqrms_L = 0;
                    freqrms_R = 0;
                    WFrms_L = 0; 
                    WFrms_R = 0;

                end
              
                track = signal_names(t);
          
                output = [output; track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, A_L, A_R, CCIR_L, CCIR_R, clicks_L, clicks_R, commonclicksa1_L, commonclicksa1_R, commonclicksb1_L, commonclicksb1_R, commonclicksa2_L, commonclicksa2_R, commonclicksb2_L, commonclicksb2_R, RMSclicks_L, RMSclicks_R, THD_L, THD_R, test_freq_L, wfreqspecamplitude_L, freqrms_L, WFrms_L, test_freq_R, wfreqspecamplitude_R, freqrms_R, WFrms_R, stereo_bleed];
               
    
            end
    end 