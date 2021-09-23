


function output = LacquerProcess(file, offset)
        [tracks, info_array] = SeperateLacquer(file, offset);
        fs = 44100;
        signal_names = tracks.keys;
        signals = tracks.values;
        output = [];
        
        
        lagdiff = info_array(1);
        normalization_L = info_array(2);
        normalization_R = info_array(3);
            for t = 1:length(signal_names)
                % for each track we need: 
                %  - RMS level
                %  - Clicks
                %  - THD
                
                % for unique tracks we need 
                %  - normalization
                %  - stereo bleed (still really unsure about this test)
                %  - wow and flutter 
                track_name = signal_names{t};
                sig = signals{t};
                csig = [];
                disp(strcat('track  ...',track_name))
    
                [csig(:,1), CLICKS_L] = ClickDetect(sig(:,1));
                [csig(:,2), CLICKS_R] = ClickDetect(sig(:,2));



               
                % [~, REFSa_L] = ClickDetect(refT(:,1));
                % [~, REFSa_R] = ClickDetect(refT(:,2));
                % [~, REFSb_L] = ClickDetect(refT(:,1));
                % [~, REFSb_R] = ClickDetect(refT(:,2));
    
                % % need to do the reference here by track 
                % commonclicksa_L = CommonClicks(CLICKS_L, REFSa_L);
                % commonclicksa_R = CommonClicks(CLICKS_R, REFSa_R);
                % commonclicksb_L = CommonClicks(CLICKS_L, REFSb_L);
                % commonclicksb_R = CommonClicks(CLICKS_R, REFSb_R);

                commonclicksa_L = 0;
                commonclicksa_R = 0;    
                commonclicksb_L = 0;
                commonclicksb_R = 0;
                lagdiff = offset*fs;
    
                
                clicks_L = length(CLICKS_L);
                clicks_R = length(CLICKS_R); 
    
                RMS_L = 20.0*log10(rms(csig(:,1)));
                RMS_R = 20.0*log10(rms(csig(:,2)));
    
                Aw = audio_Aweighting(csig(:,1));
                CCIRw = audio_CCIRweighting(csig(:,1));
    
                A_L = 20.0*log10(rms_response(Aw(1,:)));
                A_R = 20.0*log10(rms_response(Aw(2,:)));
    
                CCIR_L = 20.0*log10(avg_response(CCIRw(1,:)));
                CCIR_R = 20.0*log10(avg_response(CCIRw(2,:)));
    
                %***    DEBUG    ***%
                % figure(t); grid on;
                % plot(csig)
                % title(track_name)
                %*** DEBUG ENDS ***%
    
    
                THD_L = thd(csig(:,1),fs);
                THD_R = thd(csig(:,2),fs);
    
    
                if ismember(signal_names(t), {'1kHzL', 'sweepL', '1kHzL2', 'sweepL2'})             
                    stereo_bleed = StereoBleed(csig,1);
                elseif ismember(signal_names(t), {'1kHzR', 'sweepR', '1kHzR2', 'sweepR2'}) 
                    stereo_bleed = StereoBleed(csig,2);
                else
                    stereo_bleed = 0;
                end
    
                if ismember(signal_names(t), {'3150Hz', '3150Hz2'})
                    wow_L = WowFlutter(csig(:,1));
                    wow_R = WowFlutter(csig(:,2));
                else 
                    wow_L = 0;
                    wow_R = 0;
                end
              
                track = signal_names(t);
         
    
    
                output = [output; track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, A_L, A_R, CCIR_L, CCIR_R, clicks_L, clicks_R, commonclicksa_L, commonclicksa_R ,commonclicksb_L, commonclicksb_R, THD_L, THD_R, wow_L, wow_R, stereo_bleed];
            end
    end 