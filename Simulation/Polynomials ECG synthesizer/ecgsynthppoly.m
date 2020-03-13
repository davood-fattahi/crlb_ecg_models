function ECG=ecgsynthppoly(L,HRmean,HRdev,pct,pmean,pdev,fs,noisdev)




%%
    
%%% scaling to [0 2*pi]
beta=range(pct,'all')./(2*pi);
gamma=pct(1);
pcphase=pct./beta+gamma;

for i=1:size(pct,1)
PolyOrder=size(pmean{i,1},2)-1;    
coefs=flip(pmean{i,1});
C=zeros(PolyOrder+1);
for j=0:PolyOrder
    for k=0:j
        C(j+1,k+1)=nchoosek(j,k)*(beta^k)*(gamma^(j-k))*coefs(j+1);
    end
end
coefs=sum(C,1);
P{i,1}=flip(coefs);
end

%%%%

%%
Ts=1./fs;
HR=HRmean+HRdev*randn(1);
ECG=zeros(3,L);
phase=zeros(L,1);

for i=2:L
    omega=(HR)*2*pi + noisdev(3).*randn(1);
    phase(i)=rem(phase(i-1)+omega*Ts,2*pi)+noisdev(1)*randn(1);
    if phase(i)< phase(i-1)/2
        HR=HRmean+HRdev*randn(1);
    end
    
    k=(sum(((pcphase-phase(i))>=0),2)==1);
    p=mean(cell2mat(P(k,1)),1);

    ECG(1,i)=phase(i);
    ECG(2,i)=polyval(((p+pdev.*randn(size(pdev)))),phase(i))+noisdev(2).*randn(1);
    ECG(3,i)=omega;
end