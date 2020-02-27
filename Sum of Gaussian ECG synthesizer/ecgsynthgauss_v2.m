function ECG = ecgsynthgauss_v2(L,HRmean,HRdev,paramean,paramdev,fs,noisdev)

Ts=1./fs;
HR=HRmean+HRdev*randn(1);

%% sum of Gaussians model definition:

%%%% initial setting of variables
mdl=zeros(3,L);
mdl(:,1)=[0 ; 0 ; (HR)*2*pi];
theta=paramean+paramdev.*randn(1,15);
beta=1;
% model construction
for k=1:L
    w=1*randn(3,1).*noisdev;
    mdl(:,k+1)=processmodel_15x3_v2(mdl(:,k),Ts,beta,theta,w);
    if mdl(1,k+1)< mdl(1,k)/2
        theta=paramean+paramdev.*randn(1,15);
        HR=HRmean+HRdev*randn(1);
        mdl(3,k+1)=(HR)*2*pi;
    end
end

ECG=mdl;

