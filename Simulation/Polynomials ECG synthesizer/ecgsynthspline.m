function ECG=ecgsynthspline(L,HRmean,HRdev,p,pdev,fs,noisdev)

pmean=p.coefs;


%%
    
%%% scaling to [0 2*pi]
beta=range(p.breaks,'all')./(2*pi);
gamma=p.breaks(1);
p.breaks=p.breaks./beta+gamma;

for i=1:p.pieces
coefs=flip(pmean(i,:));
C=zeros(p.order);
for j=0:p.order-1
    for k=0:j
        C(j+1,k+1)=nchoosek(j,k)*(beta^k)*(gamma^(j-k))*coefs(j+1);
    end
end
coefs=sum(C,1);
p.coefs(i,:)=flip(coefs);
end

%%%%

%%
Ts=1./fs;
HR=HRmean+HRdev*randn(1);
ECG=zeros(3,L);
phase=zeros(L,1);
P=p.coefs;
p.coefs=polydeviate(P,1+pdev(1)*randn(1),pdev(2)*randn(1),1+pdev(3)*randn(1));
for i=2:L
    omega=(HR)*2*pi + noisdev(3).*randn(1);
    phase(i)=rem(phase(i-1)+omega*Ts,2*pi)+noisdev(1)*randn(1);
    if phase(i)< phase(i-1)/2
        HR=HRmean+HRdev*randn(1);
        p.coefs=polydeviate(P,1+pdev(1)*randn(1),pdev(2)*randn(1),1+pdev(3)*randn(1));
    end
    ECG(1,i)=phase(i);
    ECG(2,i)=ppval(p,phase(i))+(noisdev(2).*randn(1));
    ECG(3,i)=omega;
end

end








