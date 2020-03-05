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
% fvtool(d)
ecg=filtfilt(d,ecg);
figure
plot(ecg)

[~,qrs_i_raw,~]=pan_tompkin(ecg,fs,0);
peaks = false(size(ecg));
peaks(qrs_i_raw)=true;
phase=PhaseCalculation(peaks);
% [ecgM1,ECGsd,meanPhase] = MeanECGExtraction(ecg,phase,fs,0);
% [ecgM, rrM]=ecgmean(ecg,fs,'sameni',70*fs./60);
[ecgM2, rrM]=ecgmean(ecg,fs);
figure
% plot(ecgM1(170:end));
hold on
plot(ecgM2)
% k=findknots(ecgM,5,1.3);
% hold on
% plot(k,ecgM(k),'*')
p=ppolyfit(1:10000,ecgM2,3,50,25);
% for i=1:size(P,1)
%     S=polyval(P{i},tt);
%     plot(tt,S)
%     SS(i,pci)=S;
% end



