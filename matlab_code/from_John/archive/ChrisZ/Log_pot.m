%Log potentiometer study
%showing how a parabolic response (two ganged linear pots) is loglike
clear all; clc;close all;
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

try
    pkg load signal; %for Octave
catch
end
%% -----------------------------------------------------------------
x=0:0.01:1;
yp=x.^2;
a=0.45;
yexp=exp(x/a)/exp(1/a);
yexp=(yexp-yexp(1))/(1-yexp(1));

figure(1)
plot(x,yp,'b');
grid on;hold on;
plot(x,yexp,'r')
legend('pot','true log')
title('pot rotation curves')
xlabel('fraction')



disp('-----------------------------finished------------------------')


