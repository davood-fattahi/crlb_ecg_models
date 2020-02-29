clear; close all; clc;

%%%% mean parameters:
a=[2.2 -8 30  -4   2 ]; % Gaussian Amplitudes
b=[0.08 0.03  0.02  0.03  0.14]*pi; % Gaussian width
c=[-pi/3 -pi/12 0 pi/15 5*pi/6]+3*pi/4; % Gaussian centers
paramean=  [a  b  c]; % mean of parameters

%%%%% variance of parameters:
dev_a=[0.001 0.0005 0.0005 0.0005 0.001].*a;
dev_b=[0.0001 0 0 0 0.0001].*b;
dev_c=[0.0001 0 0 0 0.0001].*c;

paramdev=[dev_a dev_b dev_c]; % deviation of parameters
% paramdev=0.001*paramean;

%%%% variance of observation noise
noisdev=[0 0.1 .03]'; % deviation of observation noise


HRmean=1.2; % Heart rate per sec 
HRdev=0.1; % deviation of Heart rate
fs=1000;
ECG = ecgsynthgauss_v1(10*fs,HRmean,HRdev,paramean,paramdev,fs,noisdev);
figure('Units','normalized','OuterPosition',[0 .5 1 .35])
plot(ECG(1,:))
hold on
plot(ECG(2,:))
hold on
plot(ECG(3,:))
legend('phase','ECG','angular velocity')