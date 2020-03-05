clear
close all
clc
load '..\Signal.mat'
ecg=S(:,1);
fs=400;
t=530:850;
s=ecg(t);
s=s-mean(s);
s=s+.1*rand(size(s));

subplot(311)
plot(s)
po=2;
vc=[];
for j=0:po
    vc(j+1)=nchoosek(po,j);
end
Vc=sqrt(sum(vc.^2));  
D=diff(s,po)./Vc;
subplot(312)
plot(D)   
 subplot(313)
 plot([zeros(po,1) ;D(:)].*s);
