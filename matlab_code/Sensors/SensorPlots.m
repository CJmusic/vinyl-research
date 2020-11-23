clear all;close all;clc
set(0,'DefaultLineLineWidth',0.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

% Tbl = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137_SensorTable.csv');
% Tbl(196,:) = [];

Tbl = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000_SensorTable.csv');
Tbl(76,:) = [];

plot_sensor(Tbl, 'PressingNumber', 'maxPressForce_Ton', 'minPressForce_Ton')
plot_sensor(Tbl, 'PressingNumber', 'maxMouldSteamIn_F', 'minMouldSteamIn_F')
plot_sensor(Tbl, 'PressingNumber', 'maxMouldSteamOutTop_F', 'minMouldSteamOutTop_F')
plot_sensor(Tbl, 'PressingNumber', 'maxMouldSteamOutBottom_F', 'minMouldSteamOutBottom_F')
plot_sensor(Tbl, 'PressingNumber', 'maxExtruderFeedthroatTemp_F', 'minExtruderFeedthroatTemp_F')
plot_sensor(Tbl, 'PressingNumber', 'maxExtruderBarrelZone1Temp_F', 'minExtruderBarrelZone1Temp_F')
plot_sensor(Tbl, 'PressingNumber', 'maxExtruderBarrelZone2Temp_F', 'minExtruderBarrelZone2Temp_F')
plot_sensor(Tbl, 'PressingNumber', 'maxExtruderBarrelZone3Temp_F', 'minExtruderBarrelZone3Temp_F')
plot_sensor(Tbl, 'PressingNumber', 'maxExtruderDieZoneTemp_F', 'minExtruderDieZoneTemp_F')
plot_sensor(Tbl, 'PressingNumber', 'maxExtruderPremouldTemp_F', 'minExtruderPremouldTemp_F')
plot_sensor(Tbl, 'PressingNumber', 'maxExtruderMeltTemp_F', 'minExtruderMeltTemp_F')

plot_fftsensor(Tbl,'maxExtruderBarrelZone1Temp_F', 'Spectrum of max extruder barrel zone 1 temperature')
plot_fftsensor(Tbl,'maxExtruderBarrelZone2Temp_F', 'Spectrum of max extruder barrel zone 2 temperature')
plot_fftsensor(Tbl,'maxExtruderDieZoneTemp_F', 'Spectrum of max extruder die zone temperature')
plot_fftsensor(Tbl,'maxExtruderPremouldTemp_F', 'Spectrum of max extruder premould Temperature')
plot_fftsensor(Tbl,'maxExtruderMeltTemp_F', 'Spectrum of extruder melt  Temperature')
plot_fftsensor(Tbl,'maxMouldSteamIn_F', 'Spectrum of max steam in temperature')
plot_fftsensor(Tbl,'maxMouldSteamOutTop_F', 'Spectrum of max steam out temperature for the top mould')
plot_fftsensor(Tbl,'maxMouldSteamOutBottom_F', 'Spectrum of max steam out temperature for the bottom mould')



plot_fftsensor(Tbl,'minExtruderBarrelZone1Temp_F', 'Spectrum of min extruder barrel zone 1 temperature')
plot_fftsensor(Tbl,'minExtruderBarrelZone2Temp_F', 'Spectrum of min extruder barrel zone 2 temperature')
plot_fftsensor(Tbl,'minExtruderDieZoneTemp_F', 'Spectrum of min extruder die zone temperature')
plot_fftsensor(Tbl,'minExtruderPremouldTemp_F', 'Spectrum of min extruder premould Temperature')
plot_fftsensor(Tbl,'minExtruderMeltTemp_F', 'Spectrum of extruder melt  Temperature')
plot_fftsensor(Tbl,'minMouldSteamIn_F', 'Spectrum of min steam in temperature')
plot_fftsensor(Tbl,'minMouldSteamOutTop_F', 'Spectrum of min steam out temperature for the top mould')
plot_fftsensor(Tbl,'minMouldSteamOutBottom_F', 'Spectrum of min steam out temperature for the bottom mould')

function plot_sensor(Tbl,x,y1,y2, titlestring)
    cols = Tbl.Properties.VariableNames;
    colx = find(ismember(cols, x));
    coly1 = find(ismember(cols, y1));
    coly2 = find(ismember(cols, y2));
    X = table2array(Tbl(:,colx));
    Y1 = table2array(Tbl(:,coly1));
    Y2 = table2array(Tbl(:,coly2));

    % figure(plotnum);  

    titlestring = y1(4:end-2);
    unit = extractAfter(y1,"_");
    if strcmp(unit,'F')
        Y1 = (Y1-32)*(5/9); %convert to celcius
        Y2 = (Y2-32)*(5/9); %convert to celcius
        unit = 'C';
    end

    fig = figure('Visible', 'off');
    scatter(X,Y1,'ko');
    grid on; hold on;
    scatter(X,Y2,'kx');
    legend(y1,y2)
    ylabel(strcat(titlestring,' [', unit,']'))
    xlabel(x)
    title(strcat(titlestring))
    plotname = strcat('plots/', titlestring,'.png');
    saveas(fig, plotname);

end

function plot_fftsensor(Tbl,x, titlestring)
    cols = Tbl.Properties.VariableNames;
    colx = find(ismember(cols, x));
    data = table2array(Tbl(:,colx));
    
    % w = hann(length(X));
    % X = X.*w;
    
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
    saveas(fig, plotname);

end %audio_spectrum