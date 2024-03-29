clc
clear
close all


load '..\ECG.mat'
ecg=ECG(2,:);
fs=1000;


d = designfilt('highpassiir','StopbandFrequency',1, ...
         'PassbandFrequency',2,'PassbandRipple',0.5, ...
         'StopbandAttenuation',65,'DesignMethod','butter', ...
         'SampleRate',fs);
ecg=filtfilt(d,ecg);

figure
subplot(211)
plot(ecg)

[ecgM, rrM]=ecgmean(ecg,fs);
t=(0:1:size(ecgM(:),1)-1)./fs;
subplot(212)
plot(t,ecgM)
ecgM=ecgM+randn(size(ecgM));

hold on
plot(t,ecgM,'y')
% Breaks
breaks = [0  .2  .26  .28 .3 .31 .32 .33 .34 .35 .36 .38 .4  .5 .65 .9430];

% Fit a spline of order 5
pp = splinefit(t,ecgM,breaks,5,'p');
Y=ppval(pp,t);
hold on
plot(t,Y)

L=5000; HRmean=1.2; HRdev=0.1; p=pp; pdev=[0 0 0]; noisdev=[0 0 0];
ECG=ecgsynthspline(L,HRmean,HRdev,p,pdev,fs,noisdev);
figure
plot(ECG(2,:));



%%
% % Fit a spline of order 3 with periodic boundary conditions
% pp = splinefit(x,y,breaks,3,'p');
% 
% % Constraints: y(0) = 0, y'(0) = 1 and y(3) + y"(3) = 0
% xc = [0 0 3];
% yc = [0 1 0];
% cc = [1 0 1; 0 1 0; 0 0 1];
% con = struct('xc',xc,'yc',yc,'cc',cc);
% 
% % Fit a cubic spline with 8 pieces and constraints
% pp = splinefit(x,y,8,con);
% 
% % Fit a spline of order 6 with constraints and periodicity
% pp = splinefit(x,y,breaks,con,6,'p');






