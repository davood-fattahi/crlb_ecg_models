clear all; 
close all; 
clc;

t=-4:.1:4;
a=1; m=0; sigma=.5;
x=a.*exp(-((t-m).^2)./(2.*(sigma.^2)));
figure
plot(t,x);


PO=4;
breaks=[-4 -3 -2 -1.5 -1 -.5 0 .5 1 1.5 2 3 4];
pp = splinefit(t,x,breaks,PO);
xh=ppval(pp,t);

hold on 
plot(t,xh);

amp=1;
beta=1;
gamma=3;

p=sppdeviate(pp,beta,gamma,amp);
tt=p.breaks(1):.1:p.breaks(end);
xhh=ppval(p,tt);
hold on 
plot(tt,xhh);





