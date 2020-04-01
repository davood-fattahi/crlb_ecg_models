clc
clear all
close all
t=-10:1:10;
m=0;
b=3;
N=10;
c=[-1 +3 -15 +105 -945]
X=zeros(floor(N./2)+1,size(t(:),1));
X(1,:)=1;
for n=2:2:N
    X(floor(n./2)+1,:)=(1./factorial(n)).*(c(n/2)./(b^n)).*(t-m).^n;
end
x=sum(X,1);
plot(x)




