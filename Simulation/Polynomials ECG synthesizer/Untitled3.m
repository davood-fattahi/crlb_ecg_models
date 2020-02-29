clc
clear
close all


load '..\Signal.mat'
ecg=S(:,1);
fs=400;
t=530:850;
s=ecg(t);
s=s-mean(s);
s=s+.1*rand(size(s));

figure
plot(t,s);

maxi=90; po=3; tr=.01;  wsi=[15 200];
[P, WS, PO, NN]=polyfitAMW(t,s,tr,maxi,'PolyOrder',po,wsi);


