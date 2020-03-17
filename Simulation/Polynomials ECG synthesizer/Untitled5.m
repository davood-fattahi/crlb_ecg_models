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
t=(0:1:size(ecgM(:),1)-1)./fs;
subplot(312)
plot(t,ecgM)

[p,pcs, pct]=ppolyfit(t,ecgM,5,50,20);
[S, SS]=ppolyval(p,pct,fs);
hold on 
plot(t,S)

subplot(313)
hold on
plot(t,ecgM)
[p,pcs, pct]=ppolyfit2(t,ecgM,5,100,0,0);
p=deviate(p,.8);
[S, SS]=ppolyval(p,pct,fs);
plot(t,S)


%%
% L=8000; HRmean=1.2; HRdev=0.1; pmean=p; pdev=0.35; noisdev=[0;0;0];
% ECG=ecgsynthppoly(L,HRmean,HRdev,pct,pmean,pdev,fs,noisdev);
% 
% figure('Units','normalized','OuterPosition',[0 .5 1 .35])
% plot(ECG(1,:))
% hold on
% plot(ECG(2,:))
% hold on
% plot(ECG(3,:))
% legend('phase','ECG','angular velocity')
% 
% 
% 






