function ECG=ecgsynthgauss_v3(L,HRmean,HRdev,paramean,paramdev,fs,noisdev)


Ts=1./fs;
HR=HRmean+HRdev*randn(1);
omega=(HR)*2*pi + noisdev(3).*randn(1);
theta=paramean+paramdev.*randn(1,15);
ECG=zeros(3,L);
phase=zeros(L,1);

for i=2:L
    omega=(HR)*2*pi + noisdev(3).*randn(1);
    phase(i)=rem(phase(i-1)+omega*Ts,2*pi)+noisdev(1)*randn(1);
    if phase(i)< phase(i-1)/2
        theta=paramean+paramdev.*randn(1,15);
        HR=HRmean+HRdev*randn(1);
    end

    
    %%% p wave
    dc(1)=rem((phase(i)-theta(11)+pi),2*pi)-pi; % distance from the center
    zp=theta(1)*exp(-(dc(1)^2)./(2*theta(6)^2)); % Gaussian function (differential form) 


    %%% Q wave
    dc(2)=rem((phase(i)-theta(12)+pi),2*pi)-pi;
    zq=theta(2)*exp(-(dc(2)^2)./(2*theta(7)^2));

    %%% R wave
    dc(3)=rem((phase(i)-theta(13)+pi),2*pi)-pi;
    zr=theta(3)*exp(-(dc(3)^2)./(2*theta(8)^2));

    %%% S wave
    dc(4)=rem((phase(i)-theta(14)+pi),2*pi)-pi;
    zs=theta(4)*exp(-(dc(4)^2)./(2*theta(9)^2));

    %%% T wave
    dc(5)=rem((phase(i)-theta(15)+pi),2*pi)-pi;
    zt=theta(5)*exp(-(dc(5)^2)./(2*theta(10)^2));

    ECG(1,i)=phase(i);
    ECG(2,i)=(zp+zq+zr+zs+zt)+noisdev(2)*randn(1); % updating the synthetic ECG
    ECG(3,i)=omega;
end