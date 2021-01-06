clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')

folder = '/Volumes/AUDIOBANK/audio_files/stylusforce/'; 

normalization = '/Volumes/AUDIOBANK/audio_files/stylusforce/normalization1.0g.wav';

files = dir(fullfile(folder,'*.wav'));

tbl = cell(3,3);
Tbl = cell2table(tbl);
Tbl.Properties.VariableNames = {'stylusforce', 'left', 'right'};

SF1 = [];
SF15 = [];
SF2 = [];
SF25 = [];
SF3 = [];



for i = (1:length(files)) 
    filename = files(i).name;

    if strcmp(filename, 'normalization1.0g.wav');
            [sig, fs] = audioread(strcat(folder,filename));

        %~~~~~~~~~~~~~~~~~~~~ NORMALIZATION ~~~~~~~~~~~~~~~~~~~~%

            L = 2^16;
            seg = sig(floor(length(sig)/2) - L/2:floor(length(sig)/2) + L/2 - 1,:);
            
            %~~~~ window the data
            win = flattopwin(L);
            winseg = seg.*win;
            windowfactor = 0.2155774;% 0.2155774 for flattop window, 0.5 for hann window


            [fftseg, fftfreq] = audio_spectrum(winseg/windowfactor, fs, 1, L);
            % figure(101)
            % audio_plotspectrum(fftfreq, fftseg, 'before normalization')

            w1 = 2*707/fs; w2 = 2*1404/fs;
            [b,a] = butter(4, [w1 w2]);
            segfilt = filter(b,a,seg);
            normalization = rms_response(seg);

            normalization_L = normalization(1);
            normalization_R = normalization(2);
            
            % disp(strcat('max data before...', num2str(max(data))))
            % disp(strcat('max data before dB...', num2str(20.0*log10(max(abs(fftseg))))))
            % disp(strcat('rms data before...', num2str(rms(data))))
            % disp(strcat('rms data before dB...', num2str(20.0*log10(rms(data)))))
            % disp(strcat('normalization_L...', num2str(normalization_L)))
            % disp(strcat('normalization_R...', num2str(normalization_R)))
            
            sig(:,1)=sig(:,1)./normalization_L;
            sig(:,2)=sig(:,2)./normalization_R;

    else 
        %get stylus force from file name 

        stylusforce = str2num(filename(3:5));
        [sig, fs] = audioread(strcat(folder,filename));
        sig(:,1)=sig(:,1)./normalization_L;
        sig(:,2)=sig(:,2)./normalization_R;
        sig = audio_Aweighting(sig);



        if stylusforce == 1.0;
            SF1 = [SF1; [rms(sig(length(sig)/3:2*length(sig)/3,1)),rms(sig(length(sig)/3:2*length(sig)/3,2))]];
        end
        if stylusforce == 1.5;
            SF15 = [SF15; [rms(sig(length(sig)/3:2*length(sig)/3,1)),rms(sig(length(sig)/3:2*length(sig)/3,2))]];
        end

        if stylusforce == 2.0;
            SF2 = [SF2; [rms(sig(length(sig)/3:2*length(sig)/3,1)),rms(sig(length(sig)/3:2*length(sig)/3,2))]];
        end


        if stylusforce == 2.5;
            SF25 = [SF25; [rms(sig(length(sig)/3:2*length(sig)/3,1)),rms(sig(length(sig)/3:2*length(sig)/3,2))]];
        end


        if stylusforce == 3.0;
            SF3 = [SF3; [rms(sig(length(sig)/3:2*length(sig)/3,1)),rms(sig(length(sig)/3:2*length(sig)/3,2))]];
        end



    end

end

SF = [1.0,1.5,2.0,2.5,3.0]

SF1
SF15
SF2
SF25
SF3 

MEANS = [];
MEANS = [mean(SF1); mean(SF15); mean(SF2); mean(SF25); mean(SF3)]
% MEANS = MEANS.'

STES = [std(SF1)/sqrt(3); std(SF15)/sqrt(2); std(SF2)/sqrt(3); std(SF25)/sqrt(3); std(SF3)/sqrt(3)]

fig = figure(1);
H = plot(SF, 20*log10(MEANS(:,1)),'ko')
hold on;
H = plot(SF, 20*log10(MEANS(:,2)),'kx')
% errorbar(SF, 20*log10(STES(:,1))-20*log10(MEANS(:,1)))
% errorbar(SF, 20*log10(STES(:,2))-20*log10(MEANS(:,2)))
grid on;
xlim([0.5,3.5])
% ylim([-39,-34])
legend('Left Channel', 'Right Channel')
xlabel('Stylus Force [g]')
ylabel('RMS level [dB]')
title('RMS noise vs Stylus Force')
% saveas(fig,'SFvsRMS.png')
saveas(fig,'SFvsA.png')

