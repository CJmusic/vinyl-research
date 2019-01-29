reference = 'correlation1.wav'
filename = 'correlation2.wav'


[data_file, fs] = audioread(filename);
[data_ref, fs] = audioread(reference);


%time = linspace(0,(length(data)-1)/fs,length(data));
%rotation_speed = 33.33333;%45;
%T = 60/rotation_speed; %this is the length of one groove segment
%num_segs = (floor(length(data)/fs/T))
%n_sam = round(T*fs)
%time_seg = time(1:n_sam);

data_ref = data_ref((length(data_ref)/3):length(data_ref));
data_file = data_file((length(data_file)/3):length(data_file));



size(data_ref)
size(data_file)


%size(data_ref.'(1))
%size(data_file.'(1))


%data_ref = data_ref.';
%data_file = data_file.';
%data_ref = data_ref;%(:,1);
%data_file = data_file;%(:,1);

%s1 = data_file;
%s2 = data_ref;

t_ref = (0:length(data_ref)-1)/fs;
t_file = (0:length(data_file)-1)/fs;

disp('SIZE MATTERS')
size(data_ref)
size(data_file)


figure(2)
subplot(2,1,1)
plot(t_ref,data_ref)
title('s_1')

subplot(2,1,2)
plot(t_file,data_file)
title('s_2')
xlabel('Time (s)')

[acor,lag] = xcorr(data_file,data_ref);

[~,I] = max(abs(acor));
lagDiff = lag(I)

timeDiff = lagDiff/fs

figure(3)
plot(lag,acor)
%a3 = gca;
%a3.XTick = sort([-3000:1000:3000 lagDiff]);

%cdata_ref = data_ref(-lagDiff+1:end);
%ct_ref = (0:length(cdata_ref)-1)/fs;

cdata_file = data_file(-lagDiff+1:end);
ct_file = (0:length(cdata_file)-1)/fs;


figure(1)
subplot(2,1,1)
plot(t_ref,data_ref)
title('s_1, aligned')
xlim([0,1])
subplot(2,1,2)


plot(ct_file,cdata_file)
title('s_2')
xlabel('Time (s)')
xlim([0,1])
