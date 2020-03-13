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
subplot(311)
plot(ecg)

[ecgM, rrM]=ecgmean(ecg,fs);
t=1:size(ecgM(:),1);
subplot(312)
plot(t,ecgM)

[p,pcs, pct]=ppolyfit(t,ecgM,5,50,20);
[S, SS]=ppolyval(p,pct,1);
hold on 
plot(S)

subplot(313)
hold on
plot(t,ecgM)
[p,pcs, pct]=ppolyfit2(t,ecgM,5,30,10,10);
[S, SS]=ppolyval(p,pct,1);
plot(S)





