% Record file integrity analysis using xcorr on sweep
% Sept 2020 John Vanderkooy
%  
clear all; clc;close all;
disp('----------------start of program--------------------')
set(0,'DefaultLineLinewidth',1.5)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
try
pkg load signal %for Octave
catch
end

% Chris: you must compare a-side dut only with a-side reference

%filename='27b_declick.wav';t1=956.8037;% example of bad file
%filename='r29a-Technics-full.wav';t1=957.7999;%smeared
%filename='r30a-Technics-full.wav';t1=958.73495;%smeared
filename='r30b-Technics-full.wav';t1=952.7588;%seems clean
%filename='r30b-Technics-full-10ms.wav';t1=952.7488;%seems clean
%filename='r31a-Technics-full.wav';t1=955.48471;%smeared
%filename='r31b-Technics-full.wav';t1=958.9805;%seems good
%filename='r104b-1-declicked.wav';t1=958.4144;%clean file
%filename='minpucksize1b.wav';t1=957.352;% smeared
%filename='minpucksizem6b.wav';t1=957.5021;% smeared
%filename='042820_late_test_file.wav';t1=957.2155;%*********reference file
%filename='042820_late_test_file_declicked.wav';t1=957.2155;%off timing
%filename='040318_early_test_file.wav';t1=950.7134;% really bad file
%filename='071819-A0000B000r24a.wav';t1=956.1219;% bad file

%---------------the reference file----------------------
reffilename='r104b-1-declicked.wav';t1ref=958.4144;%known clean reference
% use Audacity to find sweep start at 137s or 655s, duration 35sec
[ref,fs]=audioread(reffilename);%t1ref comes with filename
%----------find indices for analysis-----------
refend=round(t1ref*fs);%fiducial lockout groove for reference file
dutend=round(t1*fs);%fiducial lockout groove for file under test
%-----------------indices for 1st ref sweep---------
refsweep1start=round(137*fs);% 137 then 655 sec
offset=refend-refsweep1start;% time between sweep start and end marker
dutsweep1start=dutend-offset;
%-------------------indices for 2nd ref sweep----------------
refsweep2start=round(655*fs);% 137 or 655 sec
offset=refend-refsweep2start;% time between sweep start and end marker
dutsweep2start=dutend-offset;
%----------select reference portion----------
lr=1;% left or right channel
% take 35 sec of reference sweep files
ref1=ref(refsweep1start:refsweep1start+round(35*fs)-1,lr);
ref2=ref(refsweep2start:refsweep2start+round(35*fs)-1,lr);
clear ref % keeps memory manageable
%----------select dut signal portions----------
[dut,fs]=audioread(filename);% t1 comes with filename
% take 35 sec of dut sweep files
dut1=dut(dutsweep1start:dutsweep1start+round(35*fs)-1,lr);
dut2=dut(dutsweep2start:dutsweep2start+round(35*fs)-1,lr);
clear dut
% if dut has missing samples its sweep comes earlier wrt reference
N=length(dut1)
t=(0:N-1)/fs;
%------------------
figure(10)
plot(t,ref1)
grid on
xlabel('time [s]')
ylabel('ref1 sweep')
title('reference')
%------------------
figure(15)
plot(t,dut1)
grid on
xlabel('time [s]')
ylabel('dut1 sweep')
title('dut')
%-----------------crosscorrelations----------
nlags=100000;
[C1,tlag]=xcorr(dut1,ref1,nlags);
[C2,tlag]=xcorr(dut2,ref2,nlags);
%----------plot amplitude------------
figure(20)
plot(tlag/fs,C1,'b')
grid on;hold on
plot(tlag/fs,C2,'r')
legend('1st sweep','2nd sweep')
xlabel('time [s]')
ylabel('crosscorrelation')
axis([xlim ylim])
title('very coarse comparison')
%----------plot amplitude------------
figure(30)
plot(tlag/fs,C1,'b')
grid on;hold on
plot(tlag/fs,C2,'r')
legend('1st sweep','2nd sweep')
xlabel('time [s]')
ylabel('crosscorrelation')
axis([xlim/10 ylim])
title('positive offset means dut missing samples wrt reference')
%----------plot amplitude------------
figure(40)
plot(tlag/fs,C1,'b')
grid on;hold on
plot(tlag/fs,C2,'r')
legend('1st sweep','2nd sweep')
xlabel('time [s]')
ylabel('crosscorrelation')
axis([xlim/100 ylim])
title('very fine comparison')
disp('-------------end of test------------------')




