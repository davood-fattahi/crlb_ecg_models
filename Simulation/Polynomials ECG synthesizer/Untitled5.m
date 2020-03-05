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
plot(ecg)

[ecgM2, rrM]=ecgmean(ecg,fs);
t=1:size(ecgM2(:),1);

figure
hold on
plot(t,ecgM2)
[p,pcs]=ppolyfit(t,ecgM2,3,50,0);

SS=nan(size(pcs,1),size(ecgM2,2));
for i=1:size(pcs,1)
    %%% polynomials evaluation
    pci=pcs(i,1):pcs(i,2);
    tt=t(pci); tt=tt(:);
    S=polyval(p{i,1},tt);
    plot(tt,S)
    SS(i,pci)=S;
end
SS=mean(SS,1,'omitnan');
figure
plot(t,SS)


    
    

figure
hold on
plot(t,ecgM2)
[p,pcs]=ppolyfit2(t,ecgM2,5,40,20,20);

SS=nan(size(pcs,1),size(ecgM2,2));
for i=1:size(pcs,1)
    %%% polynomials evaluation
    pci=pcs(i,1):pcs(i,2);
    tt=t(pci); tt=tt(:);
    S=polyval(p{i,1},tt);
    plot(tt,S)
    SS(i,pci)=S;
end
SS=mean(SS,1,'omitnan');
figure
plot(t,SS)




