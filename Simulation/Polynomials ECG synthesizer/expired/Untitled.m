


close all
clear
clc
po=3;
do=2;

for j=0:do
vc(j+1)=nchoosek(do,j);
end
vc=sum(vc.^2);  



tt=1:250;
s=tt.^po-tt.^(po-1)+tt.^(po-2);
ss=tt.^2;
ss=[flip(ss) s];
ss=ss+randn(size(ss));
ss(end-240:end)=[];
D=diff(ss,do)./(vc^.5);
var(D)
plot(ss);
figure
plot(D);
