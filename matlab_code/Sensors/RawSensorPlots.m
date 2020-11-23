clear all;close all;clc
set(0,'DefaultLineLineWidth',0.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')

% Tbl = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/rawpress_data/120719_SensorValues.csv');
Tbl = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/rawpress_data/121918_SensorValues.csv');



% Tbl = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/rawpress_data/120719_SensorValues.csv');
Tbl(3000:end,:) = [];



plot_rawsensor(Tbl, 'RecordTimeStamp', 'PressPosition_Inches', 'Press position')
plot_rawsensor(Tbl, 'RecordTimeStamp', 'PressForce_Ton', 'Press force')
plot_rawsensor(Tbl, 'RecordTimeStamp', 'MouldSteamIn_F', 'Mould steam in')
plot_rawsensor(Tbl, 'RecordTimeStamp', 'MouldSteamOutTop_F', 'Mould steam out top')
plot_rawsensor(Tbl, 'RecordTimeStamp', 'MouldSteamOutBottom_F', 'Mould steam out bottom')
plot_rawsensor(Tbl, 'RecordTimeStamp', 'ExtruderFeedthroatTemp_F', 'Extruder feedthroat')
plot_rawsensor(Tbl, 'RecordTimeStamp', 'ExtruderBarrelZone1Temp_F', 'minExtruderBarrelZone1Temp_F')
plot_rawsensor(Tbl, 'RecordTimeStamp', 'ExtruderBarrelZone2Temp_F', 'minExtruderBarrelZone2Temp_F')
plot_rawsensor(Tbl, 'RecordTimeStamp', 'ExtruderDieZoneTemp_F', 'minExtruderDieZoneTemp_F')
plot_rawsensor(Tbl, 'RecordTimeStamp', 'ExtruderPremouldTemp_F', 'minExtruderPremouldTemp_F')
plot_rawsensor(Tbl, 'RecordTimeStamp', 'ExtruderMeltTemp_F', 'minExtruderMeltTemp_F')

plot_fftsensor(Tbl,'ExtruderBarrelZone1Temp_F', 'Spectrum of extruder barrel zone 1 temperature')
plot_fftsensor(Tbl,'ExtruderBarrelZone2Temp_F', 'Spectrum of extruder barrel zone 2 temperature')
plot_fftsensor(Tbl,'ExtruderDieZoneTemp_F', 'Spectrum of extruder die zone temperature')
plot_fftsensor(Tbl,'ExtruderPremouldTemp_F', 'Spectrum of extruder premould temperature')
plot_fftsensor(Tbl,'ExtruderMeltTemp_F', 'Spectrum of extruder melt temperature')
plot_fftsensor(Tbl,'MouldSteamIn_F', 'Spectrum of steam in temperature')
plot_fftsensor(Tbl,'MouldSteamOutTop_F', 'Spectrum of steam out temperature for the top mould')
plot_fftsensor(Tbl,'MouldSteamOutBottom_F', 'Spectrum of steam out temperature for the bottom mould')



function plot_rawsensor(Tbl,x,y1, titlestring)
    cols = Tbl.Properties.VariableNames;
    colx = find(ismember(cols, x));
    coly1 = find(ismember(cols, y1));
    X = table2array(Tbl(:,colx));
    Y1 = table2array(Tbl(:,coly1));

    titlestring = y1(1:end-2);
    unit = extractAfter(y1,"_");
    if strcmp(unit,'F')
        Y1 = (Y1-32)*(5/9); %convert to celcius
        unit = 'C';
    end

    fig = figure('Visible', 'off');
    scatter(X,Y1,'ko');
    grid on; hold on;

    ylabel(strcat(titlestring,' [', unit,']'))
    xlabel(x)
    % title(titlestring);
    plotname = strcat('plots/', titlestring,'.png');
    saveas(fig, plotname);

end

function plot_fftsensor(Tbl,x, titlestring)
    cols = Tbl.Properties.VariableNames;
    colx = find(ismember(cols, x));
    data = table2array(Tbl(:,colx));
    
    w = hann(length(data));
    data = data.*w;
    
    fs = 25; %the sampling rate of the press sensors
    start_sam = 1;
    n_sam = length(data);
    
    
    data_fft = (fft(data(start_sam:start_sam+n_sam-1, :))/n_sam);
    data_fft = data_fft(1:floor(n_sam/2)+1,:);
    data_fft(2:end-1,:) = 2*data_fft(2:end-1,:);
    
    %fix DC
    data_fft(1,:) = 0.5*data_fft(1,:);
    % fix Nyquist
    if floor(n_sam/2)==n_sam/2
        data_fft(floor(n_sam/2+1)) = 0.5*data_fft(floor(n_sam/2+1));% fix Nyquist
    end

    freq = fs*(0:(n_sam/2))/n_sam;
    % disp('finished audio_spectrum')
    
    plotname = strcat('plots/',titlestring,'.png');
    fig = figure('Visible', 'off');
    plot(freq, abs(data_fft), 'k') 
    grid on 
    set(gca, 'XScale', 'log');
    title(titlestring)
    xlabel('Frequency (Hz)')
    ylabel('[C^-1]')
    ylim([0,10])
    saveas(fig, plotname);

end %audio_spectrum

