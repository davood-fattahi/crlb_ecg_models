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
[Ra, Ri]=max(ecgM);
t=Ri-6:Ri+6;
tt=t; tt(7)=[];
B=(ecgM(tt)-Ra);
A=((tt-Ri).^2);
b=B*A'*((A*A')^(-1));

P=[b -2*b*Ri b*Ri^2+Ra];
S=polyval(P,t);
hold on
plot(t,S)

