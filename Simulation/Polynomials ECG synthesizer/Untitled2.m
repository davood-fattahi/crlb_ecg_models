clc
clear
close all


load '..\Signal.mat'
ecg=S(:,1);
fs=400;
ecg=ecg(1:500000);

% ecgM=ecgmean(ecg,fs,'sameni',70*fs./60);
[ecgM, rrM]=ecgmean(ecg,fs);
figure
plot(ecgM);

k=ecgknots(ecgM,3,0);