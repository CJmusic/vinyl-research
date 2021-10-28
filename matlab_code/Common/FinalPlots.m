clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/Code/vinyl-research/matlab_code/')

A0000B0000 = readtable('/Users/cz/Code/vinyl-research/matlab_code/A0000B0000/A0000B0000.csv');
A0137B0137 = readtable('/Users/cz/Code/vinyl-research/matlab_code/A0137B0137/A0137B0137.csv');

A0000B0000.pressing = A0000B0000.pressing_SettingsTable;
% A0000B0000.pressing_SettingsTable = [];
A0000B0000.Properties.VariableNames([63]) = {'freqrms_R'};
A0000B0000.freqrms_R = zeros(height(A0000B0000),1);


head(A0000B0000)
head(A0137B0137)


BOTH = [A0000B0000; A0137B0137];
writetable(BOTH,'BOTH.csv');
% cols = Tbl.Properties.VariableNames;

% A0000B0000 = BOTH(strcmp(BOTH.pressing,'A0000B0000'));

% A0137B0137 = BOTH(strcmp(BOTH.pressing,'')

% % 
% % A0000B0000 = removevars(A0000B0000,{'pressing_SettingsTable'});

% T.Properties.VariableNames{58} = 'pressing';
disp('transition')
firsttransA_La = get_stats(A0000B0000, 'transition', 'a', 'A_L')
firsttransA_Ra = get_stats(A0000B0000, 'transition', 'a', 'A_R')
firsttransA_Lb = get_stats(A0000B0000, 'transition', 'b', 'A_L')
firsttransA_Rb = get_stats(A0000B0000, 'transition', 'b', 'A_R')
secondtransA_La = get_stats(A0137B0137, 'transition', 'a', 'A_L')
secondtransA_R = get_stats(A0137B0137, 'transition', 'a', 'A_R')
secondtransA_Lb = get_stats(A0137B0137, 'transition', 'b', 'A_L')
secondtransA_R = get_stats(A0137B0137, 'transition', 'b', 'A_R')

disp('quiet')
firstquietA_La = get_stats(A0000B0000, 'quiet', 'a', 'A_L')
firstquietA_Rb = get_stats(A0000B0000, 'quiet', 'a', 'A_R')
firstquietA_Lb = get_stats(A0000B0000, 'quiet', 'b', 'A_L')
firstquietA_Rb = get_stats(A0000B0000, 'quiet', 'b', 'A_R')
secondquietA_La = get_stats(A0137B0137, 'quiet', 'a', 'A_L')
secondquietA_Ra = get_stats(A0137B0137, 'quiet', 'a', 'A_R')
secondquietA_Lb = get_stats(A0137B0137, 'quiet', 'b', 'A_L')
secondquietA_Rb = get_stats(A0137B0137, 'quiet', 'b', 'A_R')

disp('quiet2')
firstquiet2A_La = get_stats(A0000B0000, 'quiet2', 'a', 'A_L')
firstquiet2A_Ra = get_stats(A0000B0000, 'quiet2', 'a', 'A_R')
firstquiet2A_Lb = get_stats(A0000B0000, 'quiet2', 'b', 'A_L')
firstquiet2A_Rb = get_stats(A0000B0000, 'quiet2', 'b', 'A_R')
secondquiet2A_La = get_stats(A0137B0137, 'quiet2', 'a', 'A_L')
secondquiet2A_Ra = get_stats(A0137B0137, 'quiet2', 'a', 'A_R')
secondquiet2A_Lb = get_stats(A0137B0137, 'quiet2', 'b', 'A_L')
secondquiet2A_Rb = get_stats(A0137B0137, 'quiet2', 'b', 'A_R')

A_LstatsTable(1,:) = struct2table(firsttransA_La)
A_LstatsTable(2,:) = struct2table(firsttransA_Ra)
A_LstatsTable(3,:) = struct2table(firsttransA_Lb)
A_LstatsTable(4,:) = struct2table(firsttransA_Rb)
A_LstatsTable(5,:) = struct2table(firsttransA_La)
A_LstatsTable(6,:) = struct2table(firsttransA_Ra)
A_LstatsTable(7,:) = struct2table(firsttransA_La)
A_LstatsTable(8,:) = struct2table(firsttransA_La)
A_LstatsTable(9,:) = struct2table(firsttransA_La)
A_LstatsTable(10,:) = struct2table(firsttransA_La)
A_LstatsTable(11,:) = struct2table(firsttransA_La)
A_LstatsTable(12,:) = struct2table(firsttransA_La)
A_LstatsTable(13,:) = struct2table(firsttransA_La)
A_LstatsTable(14,:) = struct2table(firsttransA_La)
A_LstatsTable(15,:) = struct2table(firsttransA_La)
A_LstatsTable(16,:) = struct2table(firsttransA_La)
A_LstatsTable(17,:) = struct2table(firsttransA_La)

plot_scatter2ylim(BOTH, 'quiet', 'a', 'PressingNumber', 'RMS_L', 'RMS_R','RMS vs. pressing number side a', 'BOTHPressingNumberRMSa',[-36,-24])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'PressingNumber', 'RMS_L', 'RMS_R','RMS  vs. pressing number side b', 'BOTHPressingNumberRMSb',[-36,-24])

plot_scatter2ylim(BOTH, 'quiet', 'a', 'PressingNumber', 'A_L', 'A_R','A-weighted RMS vs. pressing number side a', 'BOTHPressingNumberARMSa',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'PressingNumber', 'A_L', 'A_R','A-weighted RMS  vs. pressing number side b', 'BOTHPressingNumberARMSb',[-58,-40])


plot_scatter2ylim(BOTH, 'quiet', 'a', 'PressingNumber', 'CCIR_L', 'CCIR_R','CCIR-weighted RMS vs. pressing number side a', 'BOTHPressingNumberCCIRa',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'PressingNumber', 'CCIR_L', 'CCIR_R','CCIR-weighted RMS  vs. pressing number side b', 'BOTHPressingNumberCCIRb',[-58,-40])


plot_scatter2ylim(BOTH, '1kHz2', 'a', 'PressingNumber', 'THD_L', 'THD_R','Total harmonic distortion  vs. pressing number side a', 'BOTHPressingNumberTHDa',[-60,-35])
plot_scatter2ylim(BOTH, '1kHz2', 'b', 'PressingNumber', 'THD_L', 'THD_R','Total harmonic distortion vs. pressing number side b', 'BOTHPressingNumberTHDb', [-60,-35])
plot_scatter2ylim(BOTH, '1kHz2', 'a', 'PressingNumber', 'THD_L', 'THD_R','Total harmonic distortion  vs. pressing number side a', 'BOTHPressingNumberTHD2a',[-60,-35])
plot_scatter2ylim(BOTH, '1kHz2', 'b', 'PressingNumber', 'THD_L', 'THD_R','Total harmonic distortion vs. pressing number side b', 'BOTHPressingNumberTHD2b',[-60,-35])

plot_scatter4(BOTH, '1kHz', '1kHz2', 'a', 'PressingNumber', 'THD_L', 'THD_R', 'Total harmonic distortion side a', 'BOTHTHD1kHz2tracksa')
plot_scatter4(BOTH, '1kHz', '1kHz2', 'b', 'PressingNumber', 'THD_L', 'THD_R', 'Total harmonic distortion side b', 'BOTHTHD1kHz2tracksb')

plot_scatter2ylimtracks(BOTH, '3150Hz', '3150Hz2', 'a', 'PressingNumber', 'WFrms_L', 'Wow frequency variation side a', 'BOTHWow2tracksa',[0,10])
plot_scatter2ylimtracks(BOTH, '3150Hz', '3150Hz2', 'b', 'PressingNumber', 'WFrms_L', 'Wow frequency variation side b', 'BOTHWow2tracksb',[0,10])


plot_scatter2ylim(BOTH, 'quiet', 'a', 'PressingNumber', 'clicks_L', 'clicks_R','Number of clicks  vs. pressing number side a', 'BOTHPressingNumberCLICKSa',[0,30])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'PressingNumber', 'clicks_L', 'clicks_R','Number of clicks vs. pressing number side b', 'BOTHPressingNumberCLICKSb',[0,30])



% plot_scatter2ylim(BOTH, '3150Hz', 'a', 'PressingNumber', 'wow_L', 'wow_R','Wow  vs. pressing number side a', 'BOTHPressingNumberWOWa')
% plot_scatter2ylim(BOTH, '3150Hz', 'b', 'PressingNumber', 'wow_L', 'wow_R','Wow  vs. pressing number side b', 'BOTHPressingNumberWOWb')
% plot_scatter(BOTH, '3150Hz', 'a', 'PressingNumber', 'centreholeoffset','Centre hole offset vs. pressing number side a', 'BOTHPressingNumberCHOa')
% plot_scatter(BOTH, '3150Hz', 'b', 'PressingNumber', 'centreholeoffset','Centre hole offset vs. pressing number side b', 'BOTHPressingNumberCHOb')


plot_scatter(BOTH, '1kHzL', 'a', 'PressingNumber', 'stereo_bleed','Stereo bleed vs pressing number side a', 'BOTHPressingNumberSBa')
plot_scatter(BOTH, '1kHzL', 'b', 'PressingNumber', 'stereo_bleed','Stereo bleed vs pressing number side b', 'BOTHPressingNumberSBb')


plot_scatter2ylim(BOTH, 'quiet', 'a', 'maxMouldSteamIn_F', 'RMS_L', 'RMS_R','RMS vs. maximum mould steam in temperature side a', 'BOTHmaxMouldSteamInRMSa',[-36,-24])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'maxMouldSteamIn_F', 'RMS_L', 'RMS_R','RMS vs. maximum mould steam in temperature side b', 'BOTHmaxMouldSteamINRMSb',[-36,-24])
plot_scatter2ylim(BOTH, 'quiet', 'a', 'maxMouldSteamIn_F', 'A_L', 'A_R','A-weighted RMS vs. maximum mould steam in temperature side a', 'BOTHmaxMouldSteamInARMSa',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'maxMouldSteamIn_F', 'A_L', 'A_R','A-weighted RMS vs. maximum mould steam in temperature side b', 'BOTHmaxMouldSteamINARMSb',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'a', 'maxMouldSteamIn_F', 'CCIR_L', 'CCIR_R','A-weighted RMS vs. maximum mould steam in temperature side a', 'BOTHmaxMouldSteamInCCIRRMSa',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'maxMouldSteamIn_F', 'CCIR_L', 'CCIR_R','A-weighted RMS vs. maximum mould steam in temperature side b', 'BOTHmaxMouldSteamINCCIRRMSb',[-58,-40])

plot_scatter2ylim(BOTH, 'quiet', 'a', 'maxMouldSteamIn_F', 'clicks_L', 'clicks_R','Number of clicks vs. maximum mould steam in temperature side a', 'BOTHmaxMouldSteamInCLICKSa',[0,30])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'maxMouldSteamIn_F', 'clicks_L', 'clicks_R','Number of clicks vs. maximum mould steam in temperature side b', 'BOTHmaxMouldSteamInCLICKSb',[0,30])


plot_scatter2ylim(BOTH, 'quiet', 'a', 'maxExtruderPremouldTemp_F', 'A_L', 'A_R','A-weighted RMS vs. maximum extruder premould temperature  side a', 'BOTHmaxExtruderPremouldTempARMSa',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'maxExtruderPremouldTemp_F', 'A_L', 'A_R','A-weighted RMS vs. maximum extruder premould temperature side b', 'BOTHmaxmaxExtruderPremouldTempARMSb',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'a', 'maxExtruderPremouldTemp_F', 'CCIR_L', 'CCIR_R','CCIR-weighted RMS vs. maximum extruder premould temperature  side a', 'BOTHmaxExtruderPremouldTempCCIRRMSa',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'maxExtruderPremouldTemp_F', 'CCIR_L', 'CCIR_R','CCIR-weighted RMS vs. maximum extruder premould temperature side b', 'BOTHmaxmaxExtruderPremouldTempCCIRRMSb',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'a', 'maxExtruderPremouldTemp_F', 'RMS_L', 'RMS_R','RMS vs. maximum extruder premould temperature  side a', 'BOTHmaxExtruderPremouldTempRMSa',[-36,-24])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'maxExtruderPremouldTemp_F', 'RMS_L', 'RMS_R','RMS vs. maximum extruder premould temperature  side b', 'BOTHmaxmaxExtruderPremouldTempRMSb',[-36,-24])
plot_scatter2ylim(BOTH, 'quiet', 'a', 'maxExtruderPremouldTemp_F', 'clicks_L', 'clicks_R','Number of clicks vs. maximum extruder premould temperature side a', 'BOTHmaxExtruderPremouldTempCLICKSa',[0,30])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'maxMouldSteamIn_F', 'clicks_L', 'clicks_R','Number of clicks vs. maximum extruder premould temperature side b', 'BOTHmaxmaxExtruderPremouldTempCLICKSb',[0,30])

plot_scatter2ylim(BOTH, 'quiet', 'a', 'minMouldSteamOutBottom_F', 'A_L', 'A_R','Minimum steam out temperature of the bottom mould vs. A-weighted RMS side a', 'BOTHminMouldSteamOutBottomARMSa',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'minMouldSteamOutBottom_F', 'A_L', 'A_R','Minimum steam out temperature of the bottom mould vs. A-weighted RMS side b', 'BOTHminMouldSteamOutBottomARMSb',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'a', 'minMouldSteamOutBottom_F', 'CCIR_L', 'CCIR_R','Minimum steam out temperature of the bottom mould vs. CCIR-weighted RMS side a', 'BOTHminMouldSteamOutBottomARMSa',[-58,-40])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'minMouldSteamOutBottom_F', 'CCIR_L', 'CCIR_R','Minimum steam out temperature of the bottom mould vs. CCIR-weighted RMS side b', 'BOTHminMouldSteamOutBottomCCIRRMSb',[-58,-40])

plot_scatter2ylim(BOTH, 'quiet', 'a', 'minMouldSteamOutTop_F', 'RMS_L', 'RMS_R','Minimum steam out temperature of the top mould vs. RMS side a', 'BOTHminMouldSteamOutTopRMSa',[-36,-24])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'minMouldSteamOutBottom_F', 'RMS_L', 'RMS_R','Minimum steam out temperature of the top mould vs. number of clicks side b', 'BOTHminMouldSteamOutTopRMSb',[-36,-24])
plot_scatter2ylim(BOTH, 'quiet', 'a', 'minMouldSteamOutTop_F', 'clicks_L', 'clicks_R','Minimum steam out temperature of the top mould vs. number of clicks side a', 'BOTHminMouldSteamOutToCLICKSa',[0,30])
plot_scatter2ylim(BOTH, 'quiet', 'b', 'minMouldSteamOutBottom_F', 'clicks_L', 'clicks_R','Minimum steam out temperature of the top mould vs. number of clicks side b', 'BOTHminMouldSteamOutTopCLICKSb',[0,30])



plot_scatter2ylim(BOTH, '1kHz', 'a', 'PressingNumber', 'RMS_L', 'RMS_R','Pressing number vs. A-weighted RMS', 'BOTH1kHzRMSa', [-36,0])
plot_scatter2ylim(BOTH, '1kHz', 'b', 'PressingNumber', 'RMS_L', 'RMS_R','RMS level of reference tone', 'BOTH1kHzRMSb', [-36,0])
plot_scatter2ylim(BOTH, '1kHz', 'a', 'PressingNumber', 'A_L', 'A_R','A-weighted RMS level of reference tone', 'BOTH1kHzARMSa', [-36,0])
plot_scatter2ylim(BOTH, '1kHz', 'b', 'PressingNumber', 'A_L', 'A_R','A weighted RMS level of reference tone', 'BOTH1kHzARMSb', [-36,0])


plot_scatter2(BOTH, '1kHz', 'a', 'maxMouldSteamOutTop_F', 'THD_L', 'THD_R','Max Mould Steam Out Top vs. THD Left Channel', 'BOTH1kHzTHDLa')
plot_scatter2(BOTH, '1kHz', 'b', 'maxMouldSteamOutTop_F', 'THD_L', 'THD_R','Max Mould Steam Out Top vs. THD Left Channel', 'BOTH1kHzTHDLb')
plot_scatter2(BOTH, '1kHz2', 'a', 'maxMouldSteamOutTop_F', 'THD_L', 'THD_R','Max Mould Steam Out Top vs. THD Left Channel', 'BOTH1kHz2THDLa')
plot_scatter2(BOTH, '1kHz2', 'b', 'maxMouldSteamOutTop_F', 'THD_L', 'THD_R','Max Mould Steam Out Top vs. THD Left Channel', 'BOTH1kHz2THDLb')




plot_scatter2ylim(BOTH, 'transition', 'a', 'PressingNumber', 'clicks_L', 'clicks_R','Number of clicks in the transition track', 'BOTHclickstransitiona',[0,30])
plot_scatter2ylim(BOTH, 'transition', 'b', 'PressingNumber', 'clicks_L', 'clicks_R','Number of clicks in the transition track', 'BOTHclickstransitionb',[0,30])



plot_scatter(BOTH, '1kHzL', 'a', 'PressingNumber', 'stereo_bleed','Stereo bleed in the 1kHzL track side a', 'BOTH1kHzStereobleeda')
plot_scatter(BOTH, '1kHzL', 'b', 'PressingNumber', 'stereo_bleed','Stereo bleed in the 1kHzL track side b', 'BOTH1kHzStereobleedb')

plot_stereohistogramstats(A0000B0000, 'quiet', 'a', 'A_L', 'A_R', 'A-weighted noise in the quiet track of the first pressing side a', 'A0000B0000quietAhista')
plot_stereohistogramstats(A0000B0000, 'quiet', 'b', 'A_L', 'A_R', 'A-weighted noise in the quiet track of the first pressing side b', 'A0000B0000quietAhistb')
plot_stereohistogramstats(A0137B0137, 'quiet', 'a', 'A_L', 'A_R', 'A-weighted noise in the quiet track of the second pressing side a', 'A0137B0137quietAhista')
plot_stereohistogramstats(A0137B0137, 'quiet', 'b', 'A_L', 'A_R', 'A-weighted noise in the quiet track of the second pressing side b', 'A0137B0137quietAhistb')

plot_stereohistogramstats(BOTH, 'quiet', 'a', 'A_L', 'A_R', 'A-weighted noise in the quiet track side a', 'BOTHquietAhista')
plot_stereohistogramstats(BOTH, 'quiet', 'b', 'A_L', 'A_R', 'A-weighted noise in the quiet track side b', 'BOTHquietAhistb')


plot_stereohistogram(BOTH, 'transition', 'a', 'A_L', 'A_R', 'A-weighted noise in the transition track side a', 'transitionAhista')
plot_stereohistogram(BOTH, 'transition', 'b', 'A_L', 'A_R', 'A-weighted noise in the transition track side b', 'transitionAhistb')
plot_barchart2(BOTH, 'transition', 'a', 'A_L', 'A_R','A weighted noise in the transition track side a', '[dB]', 'transitionpressingAa')
% plot_barchart2(BOTH, 'transition', 'b', 'A_L', 'A_R','A weighted noise in the transition track side b', '[dB]', 'transitionpressingAb')

% plot_scatter2ylim(BOTH, '3150Hz', 'a', 'PressingNumber','wow_L', 'wow_R','Peak to peak wow per record', 'wowa3150')
% plot_scatter2ylim(BOTH, '3150Hz', 'b', 'PressingNumber', 'wow_L', 'wow_R','Peak to peak wow per record', 'wowb3150')
% plot_scatter2ylim(BOTH, '3150Hz', 'a', 'PressingNumber','wow_L', 'wow_R','Peak to peak wow per record', 'wowa3150')
% plot_scatter2ylim(BOTH, '3150Hz', 'b', 'PressingNumber', 'wow_L', 'wow_R','Peak to peak wow per record', 'wowb3150')


A_L = BOTH{:,'A_L'}./20;
A_R = BOTH{:,'A_R'}./20;
A_Labs = 10.^(A_L);
A_Rabs = 10.^(A_R);

BOTH(:,'A_Labs') = num2cell(A_Labs);
BOTH(:,'A_Rabs') = num2cell(A_Rabs);


% BOTH(:,'A_Labs') = 10.^((BOTH{:,'A_L'}./20));
% BOTH(:,'A_Rabs') = 10.^((BOTH{:,'A_R'}./20));

% plot_scatter2ylim(BOTH, 'quiet', 'a', 'PressingNumber', 'A_Labs', 'A_Rabs','A-weighted RMS levels per record', 'quietPressingNumberAabsa')
% plot_scatter2ylim(BOTH, 'quiet', 'b', 'PressingNumber', 'A_Labs', 'A_Rabs','A-weighted RMS levels per record', 'quietPressingNumberAabsb')

% plot_scatter2ylim(A0000B0000, 'quiet', 'a', 'PressingNumber', 'A_L', 'A_R','A-weighted RMS levels per record', 'A0000B0000quietPressingNumberARMSa')
% plot_scatter2ylim(A0000B0000, 'quiet', 'b', 'PressingNumber', 'A_L', 'A_R','A-weighted RMS levels per record', 'A0000B0000quietPressingNumberARMSb')

function plot_clicks(Tbl,side, titlestring)
    clicktracks = {{'100Hz'     }
    {'100Hz2'    }
    {'10kHz'     }
    {'10kHz2'    }
    {'1kHz'      }
    {'1kHz2'     }
    {'1kHzL'     }
    {'1kHzL2'    }
    {'1kHzR'     }
    {'1kHzR2'    }
    {'1kHzV'     }
    {'1kHzV2'    }
    {'3150Hz'    }
    {'3150Hz2'   }
    {'quiet'     }
    {'quiet2'    }
    {'sweep'     }
    {'sweep2'    }
    {'sweepL'    }
    {'sweepL2'   }
    {'sweepR'    }
    {'sweepR2'   }
    {'sweepV'    }
    {'sweepV2'   }
    {'transition'}};
    Tbl = Tbl(strcmp(Tbl.side,side),:);

    cols = Tbl.Properties.VariableNames;
    % clicks_L = cell2table(cell(0,length(cols)));
    % clicks_R = cell2table(cell(0,length(cols)));
    clicks_L = [];
    clicks_R = [];
    % for i = (1:length(clicktracks))
    Tbl(strcmp(Tbl.track,'leadout'),:) = [];
    Tbl(strcmp(Tbl.track,'leadin'),:) = [];

    for i = (1:137)
        tbl = Tbl((Tbl.PressingNumber==i),:);
        clicks_l = sum(tbl.clicks_L);
        clicks_r = sum(tbl.clicks_R);
        % for j = 1:length(clicktracks)
        %     tbl = Tbl(strcmp(Tbl.track,clicktracks{j}),:);
        %     clicks_l = tbl.clicks_L
        %     clicks_r = tbl.clicks_R
        % end
    
        clicks_L = [clicks_L; (clicks_l)];
        clicks_R = [clicks_R; (clicks_r)];
    
    end 
    pressing_number = tbl.PressingNumber;
    size(pressing_number)
    size(clicks_L)
    size(clicks_R)
    pressing_number = (1:137);
    fig = figure('Visible', 'off');
    scatter(pressing_number,clicks_L,'ko');
    grid on; hold on;
    scatter(pressing_number,clicks_R,'kx');
    legend({'left channel', 'right channel'});
    title(titlestring);
    % xlim([0,100])
    % ylim([1,160])
    xlabel('record number')
    ylabel('number of clicks')
    plotname = strcat('plots/', titlestring,'.png');
    saveas(fig, plotname);

end



function stats = get_stats(Tbl, trackname, side, column)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, column));
    X = table2array(Tbl(:,colx));

    stats = datastats(X);

end

function plot_barchart2(Tbl, trackname, side, x1, x2,titlestring, varname, filename)
    cols = Tbl.Properties.VariableNames;
    pressruns = unique(Tbl.pressing);

    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx1 = find(ismember(cols, x1));
    colx2 = find(ismember(cols, x2));
    % X = table2array(Tbl(:,colx));

    X1 = [];
    X2 = [];
    for i = 1:length(pressruns)
        tbl = Tbl(strcmp(Tbl.pressing,pressruns{i}),:);
        x1 = mean(table2array(tbl(:,colx1)));
        x2 = mean(table2array(tbl(:,colx2)));
        X1 = [X1, x1];
        X2 = [X2, x2];
    end


    % figure(plotnum); grid on; hold on;
    fig = figure('Visible', 'off'); grid on; hold on;

    H = bar([X1, X2], 'LineWidth', 2);

    hold on;
    % H = bar(X1, 'LineWidth', 2)
    H(1).FaceColor = [0.6 0.6 0.6];
    % H(2).FaceColor = [.9 .9 .9];
    set(gca,'xticklabel',pressruns)
    ax=gca;
    ax.FontSize=8;
    ax.XTick = (1:length(pressruns))   %THIS WAY, YOU SET HOW MANY XTICKS YOU WANT FOR YOUR XTICKLABELS
    xtickangle(45)
    ylabel(varname)
    title(titlestring)
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);
end

function plot_barchart(Tbl, trackname, side, x, titlestring, varname, filename)
    cols = Tbl.Properties.VariableNames;
    pressruns = unique(Tbl.pressing);

    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, x));
    % X = table2array(Tbl(:,colx));

    X = [];
    for i = 1:length(pressruns)
        tbl = Tbl(strcmp(Tbl.pressing,pressruns{i}),:);
        x = mean(table2array(tbl(:,colx)));
        X = [X, x];
    end


    % figure(plotnum); grid on; hold on;
    fig = figure('Visible', 'off'); grid on; hold on;

    H = bar(X, 'LineWidth', 2)
    H(1).FaceColor = [0.6 0.6 0.6];
    % H(2).FaceColor = [.9 .9 .9];
    set(gca,'xticklabel',pressruns)
    ax=gca;
    ax.FontSize=8;
    ax.XTick = (1:length(pressruns))   %THIS WAY, YOU SET HOW MANY XTICKS YOU WANT FOR YOUR XTICKLABELS
    xtickangle(45)
    ylabel(varname)
    title('RMS noise in quiet track')
    saveas(figure(plotnum),'RMSquiet.png')
end

function plot_histogram(Tbl, trackname, side, x, titlestring, varname, binlims, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, x));
    X = table2array(Tbl(:,colx));

    % figure(plotnum); grid on; hold on;
    fig = figure('Visible', 'off'); grid on; hold on;

    histogram(X,50,'BinLimits',binlims)

    ylabel('number of records')
    varname
    xlabel(varname)
    title(titlestring)
    plotname = strcat('plots/', filename,'.png');
    saveas(figure(plotnum), plotname)
end


function plot_scatter(Tbl, trackname, side, x, y,titlestring, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, x));
    coly = find(ismember(cols, y));
    X = table2array(Tbl(:,colx));
    Y = table2array(Tbl(:,coly));

    % figure(plotnum);  
    fig = figure('Visible', 'off');
    scatter(X,Y,'ko');
    grid on; 
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end

function plot_scatter2(Tbl, trackname, side, x, y1, y2, titlestring, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, x));
    coly1 = find(ismember(cols, y1));
    coly2 = find(ismember(cols, y2));
    X = table2array(Tbl(:,colx));
    Y1 = table2array(Tbl(:,coly1));
    Y2 = table2array(Tbl(:,coly2));

    % figure(plotnum);  
    fig = figure('Visible', 'off');
    scatter(X,Y1,'ko');
    grid on; hold on;
    scatter(X,Y2,'kx');
    legend({'left channel', 'right channel'});
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end

function plot_scatter2ylim(Tbl, trackname, side, x, y1, y2, titlestring, filename, ylims)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx = find(ismember(cols, x));
    coly1 = find(ismember(cols, y1));
    coly2 = find(ismember(cols, y2));
    X = table2array(Tbl(:,colx));
    Y1 = table2array(Tbl(:,coly1));
    Y2 = table2array(Tbl(:,coly2));

    % figure(plotnum);  
    fig = figure('Visible', 'off');
    scatter(X,Y1,'ko');
    grid on; hold on;
    scatter(X,Y2,'kx');
    legend({'left channel', 'right channel'});
    title(titlestring);
    ylim(ylims)
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end


function plot_scatter2ylimtracks(Tbl, trackname1, trackname2, side, x, y1, titlestring, filename, YLIM)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    Tbl1 = Tbl(strcmp(Tbl.track,trackname1),:);
    Tbl2 = Tbl(strcmp(Tbl.track,trackname2),:);
    % Tbl3 = Tbl(strcmp(Tbl.track,trackname3),:);
    colx = find(ismember(cols, x));
    coly1 = find(ismember(cols, y1));
    % coly2 = find(ismember(cols, y2));
    X1 = table2array(Tbl1(:,colx));
    Y1L = table2array(Tbl1(:,coly1));
    % Y1R = table2array(Tbl1(:,coly2));
    X2 = table2array(Tbl2(:,colx));
    Y2L = table2array(Tbl2(:,coly1));
    % Y2R = table2array(Tbl2(:,coly2));

    % figure(plotnum);  
    fig = figure('Visible', 'off');
    scatter(X1,Y1L,'ko');
    grid on; hold on;
    % scatter(X1,Y1R,'kx');
    scatter(X2,Y2L,'bo');
    grid on; hold on;
    ylim(YLIM);
    % scatter(X2,Y2R,'bx');
    legend(trackname1, trackname2)
    % legend({strcat(trackname1, ' left'), strcat(trackname1, ' right'), strcat(trackname2, ' left'), strcat(trackname2, ' right')});
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end

function plot_scatter4(Tbl, trackname1, trackname2, side, x, y1, y2, titlestring, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    Tbl1 = Tbl(strcmp(Tbl.track,trackname1),:);
    Tbl2 = Tbl(strcmp(Tbl.track,trackname2),:);
    % Tbl3 = Tbl(strcmp(Tbl.track,trackname3),:);
    colx = find(ismember(cols, x));
    coly1 = find(ismember(cols, y1));
    coly2 = find(ismember(cols, y2));
    X1 = table2array(Tbl1(:,colx));
    Y1L = table2array(Tbl1(:,coly1));
    Y1R = table2array(Tbl1(:,coly2));
    X2 = table2array(Tbl2(:,colx));
    Y2L = table2array(Tbl2(:,coly1));
    Y2R = table2array(Tbl2(:,coly2));

    % figure(plotnum);  
    fig = figure('Visible', 'off');
    scatter(X1,Y1L,'ko');
    grid on; hold on;
    scatter(X1,Y1R,'kx');
    scatter(X2,Y2L,'bo');
    grid on; hold on;
    scatter(X2,Y2R,'bx');
    legend({strcat(trackname1, ' left'), strcat(trackname1, ' right'), strcat(trackname2, ' left'), strcat(trackname2, ' right')});
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end


function plot_scatterraw(X, Y, titlestring, filename)
    fig = figure('Visible', 'off');
    scatter(X,Y,'ko');
    grid on; 
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end

function plot_scatterraw2(X, Y, titlestring, filename)
    fig = figure('Visible', 'off');
    scatter(X,Y,'ko');
    grid on; 
    title(titlestring);
    plotname = strcat('plots/', filename,'.png');
    saveas(fig, plotname);

end



function plot_stereohistogramstats(Tbl, trackname, side, x1, x2, titlestring, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx1 = find(ismember(cols, x1));
    colx2 = find(ismember(cols, x2));
    data_L = table2array(Tbl(:,colx1));
    data_R = table2array(Tbl(:,colx2));


    statsL = datastats(data_L);
    statsR = datastats(data_R);    
    
    % lower_binL = statsL.mean - 10;
    % lower_binR = statsR.mean - 10;
    % upper_binL = statsL.mean + 10;
    % upper_binR = statsR.mean + 10;


    fig = figure('Visible', 'off');
    % histogram(data_L, 50,'BinLimits',[lower_binL,upper_binL], 'facecolor',[0.6 0.6 0.6])
    histogram(data_L, 50, 'facecolor',[0.6 0.6 0.6])
    hold on; grid on;
    histogram(data_R,50, 'facecolor',[0.3 0.3 0.3])
    % histogram(data_R,50,'BinLimits',[lower_binL,upper_binL], 'facecolor',[0.3 0.3 0.3])
    title(titlestring)
    legend('left channel', 'right channel')
    dim = [0.2 0.5 0.3 0.3];
    str = {strcat('left mean :',num2str(statsL.mean)),strcat('left std :',num2str(statsL.std)),strcat('right mean :',num2str(statsR.mean)),strcat('right std :',num2str(statsR.std))};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor', 'white');
    plotname = strcat('plots/',filename,'.png');
    saveas(fig, plotname);
end


function plot_stereohistogram(Tbl, trackname, side, x1, x2, titlestring, filename)
    cols = Tbl.Properties.VariableNames;
    Tbl = Tbl(strcmp(Tbl.track,trackname),:);
    Tbl = Tbl(strcmp(Tbl.side,side),:);
    colx1 = find(ismember(cols, x1));
    colx2 = find(ismember(cols, x2));
    data_L = table2array(Tbl(:,colx1));
    data_R = table2array(Tbl(:,colx2));


    statsL = datastats(data_L);
    statsR = datastats(data_R);    
    
    % lower_binL = statsL.mean - 10;
    % lower_binR = statsR.mean - 10;
    % upper_binL = statsL.mean + 10;
    % upper_binR = statsR.mean + 10;


    fig = figure('Visible', 'off');
    % histogram(data_L, 50,'BinLimits',[lower_binL,upper_binL], 'facecolor',[0.6 0.6 0.6])
    histogram(data_L, 50, 'facecolor',[0.6 0.6 0.6])
    hold on; grid on;
    histogram(data_R,50, 'facecolor',[0.3 0.3 0.3])
    % histogram(data_R,50,'BinLimits',[lower_binL,upper_binL], 'facecolor',[0.3 0.3 0.3])
    title(titlestring)
    legend('left channel', 'right channel')
    % dim = [0.2 0.5 0.3 0.3];
    % str = {strcat('left mean :',num2str(statsL.mean)),strcat('left std :',num2str(statsL.std)),strcat('right mean :',num2str(statsR.mean)),strcat('right std :',num2str(statsR.std))};
    % annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor', 'white');
    plotname = strcat('plots/',filename,'.png');
    saveas(fig, plotname);
end