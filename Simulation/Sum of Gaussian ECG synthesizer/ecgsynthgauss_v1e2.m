function ECG = ecgsynthgauss_v1e2(L,HRmean,HRdev,paramean,paramdev,fs,noisdev)

Ts=1./fs;
HR=HRmean+HRdev*randn(1);

%% sum of Gaussians model definition:

%%%% initial setting of variables
mdl=zeros(3,L);
mdl(:,1)=[0 ; 0 ; (HR)*2*pi];
theta=paramean+paramdev.*randn(1,15);
beta=1;
% model construction
for k=2:L
    mdl(1,k)=mod(mdl(1,k-1)+mdl(3,k-1)*Ts,2*pi)+noisdev(1).*randn(1);
    mdl(2,k)=processmodel_15x3_v1e2(mdl(:,k-1),Ts,beta,theta);
    if mdl(1,k)< mdl(1,k-1)/2
        theta=paramean+paramdev.*randn(1,15);
        HR=HRmean+HRdev*randn(1);
    end
    mdl(3,k)=(HR)*2*pi+noisdev(3)*randn(1);
end

ECG=mdl;
ECG(2,:)=ECG(2,:)+noisdev(2).*randn(1,L);

