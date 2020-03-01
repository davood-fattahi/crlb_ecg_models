clear all
close all
clc


range0=[1:300];

beta=2/(300-1);
gamma=-1;
range1=range0.*beta+gamma;


coefs=[1 2 1];
plot(range1,polyval(flip(coefs),range1));


%%% range denormalization
order=2;

    
for i=0:order
    for j=0:i
        C(i+1,j+1)=nchoosek(i,j)*(beta^j)*(gamma^(i-j))*coefs(i+1);
    end
end
coefs=sum(C,1);

figure
plot(range0,polyval(flip(coefs),range0));












