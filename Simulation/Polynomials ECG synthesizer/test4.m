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
ecgM=ecgM+randn(size(ecgM));
t=1:size(ecgM(:),1);
subplot(312)
plot(t,ecgM)

afp =ecgfidapprox(ecgM,fs,72./60,'widest'); 
afp(3,5)=afp(3,6); afp(:,6)=[];
hold on
plot(afp(2,:),ecgM(afp(2,:)),'r*')
plot(afp(1,:),ecgM(afp(1,:)),'ro')
plot(afp(3,:),ecgM(afp(3,:)),'ro')

afp(2,:)=[];
[p,pcs, pct]=ppolyfit2(t,ecgM,[3 7 6 7 3],afp',fs./200);
[S, SS]=ppolyval(p,pct,1);
subplot(313)
plot(t,ecgM);
hold on 
plot(S)