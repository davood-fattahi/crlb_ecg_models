clear 
close all
clc

NumGauss=20;
t=.004:.004:2;
SNRdB=(-20:5:60); % vector
SNR=10.^(SNRdB/10);
pMean=[1 .1 1];
pStd=[.1 .01 .02];

PrParams=zeros(3,NumGauss);
g=zeros(length(t),NumGauss);
for j=1:NumGauss
    PrParams(:,j)=pMean+pStd.*randn(1,3);
    g(:,j)=gauss(t,PrParams(:,j));
end

PrMean=mean(PrParams,2);
PrCov=cov(PrParams');
VarSig=var(g(:));

MlParams =zeros(3,NumGauss, length(SNR));
MlEr =zeros(size(MlParams ));     

BysParams =zeros(3,NumGauss, length(SNR));
BysEr =zeros(size(BysParams ));     

VarNoiseAdd =zeros(size(SNR)); 

MlCov =zeros(3,3,length(SNR));
BysCov =zeros(3,3,length(SNR));

DetMlCov =zeros(size(SNR)); 
DetBysCov =zeros(size(SNR)); 

h = waitbar(0,'Estimating parameters for noisy signals, please wait ...');

for k=1:length(SNR)
    waitbar(k/length(SNR))
        VarNoiseAdd(k)=VarSig /SNR(k); 
        signalN =g+randn(size(g)).*sqrt(VarNoiseAdd(k));
        for j=1:NumGauss
            p0=[.8 .2 .7];
            lb=[-inf 0 t(1)]; ub=[+inf range(t) t(end)];
            %%% for ML, fit gaussians on each noisy segment
            MlParams(:,j,k)=GausFit(t,signalN(:,j) ,p0,lb(:),ub(:));
            %%% for Bayesian, fit gaussians on each noisy segment
            options = struct('Display','final-detailed' ,'MaxFunctionEvaluations',10000, 'MaxIterations', 10000, 'ConstraintTolerance', 1e-20, 'OptimalityTolerance', 1e-20, 'StepTolerance', 1e-20);
            BysParams (:,j,k)=GausFit(t,signalN(:,j), p0, lb(:), ub(:), PrMean ,PrCov, VarNoiseAdd(k),options);
        end
    MlEr(:,:,k)=MlParams (:,:,k)-PrParams;
    BysEr(:,:,k)=BysParams (:,:,k)-PrParams;
    MlCov(:,:,k)=cov(MlEr(:,:,k)');
    DetMlCov(k)=det(MlCov (:,:,k));
    BysCov(:,:,k)=cov(BysEr(:,:,k)');
    DetBysCov(k)=det(BysCov (:,:,k));
end
close(h)

MlCrlbNum=zeros(3, 3, NumGauss, length(SNRdB));
MlFINum=zeros(3, 3, NumGauss, length(SNRdB));
BysCrlbNum=zeros(3, 3, length(SNRdB));
MlCrlbNum_AvrgdFI=zeros(3, 3, length(SNRdB));
DetMlCrlbNum=zeros(length(SNRdB),1);
DetBysCrlbNum=zeros(length(SNRdB),1);

MlCrlbAppr=zeros(3, 3, NumGauss, length(SNRdB));
MlFIAppr=zeros(3, 3, NumGauss, length(SNRdB));
BysCrlbAppr=zeros(3, 3, length(SNRdB));
MlCrlbAppr_AvrgdFI=zeros(3, 3, length(SNRdB));
DetMlCrlbAppr=zeros(length(SNRdB),1);
DetBysCrlbAppr=zeros(length(SNRdB),1);

for k=1:length(SNRdB)
    for j=1:NumGauss
       [MlCrlbNum(:,:,j,k), MlFINum(:,:,j,k)]=GaussCrlbNumeric(t, PrParams(:,j), VarNoiseAdd(k), 1 );
       [MlCrlbAppr(:,:,j,k), MlFIAppr(:,:,j,k)]=GaussCRLB(PrParams(1,j), PrParams(2,j), 250, VarNoiseAdd(k) );   
    end
    MlCrlbNum_AvrgdFI(:,:,k)=inv(mean(MlFINum(:,:,:,k),3));
    BysCrlbNum(:,:,k)=inv(mean(MlFINum(:,:,:,k),3)+inv(PrCov));

    DetMlCrlbNum(k)=det(MlCrlbNum_AvrgdFI(:,:,k));
    DetBysCrlbNum(k)=det(BysCrlbNum(:,:,k));
    
    MlCrlbAppr_AvrgdFI(:,:,k)=inv(mean(MlFIAppr(:,:,:,k),3));
    BysCrlbAppr(:,:,k)=inv(mean(MlFIAppr(:,:,:,k),3)+inv(PrCov));

    DetMlCrlbAppr(k)=det(MlCrlbAppr_AvrgdFI(:,:,k));
    DetBysCrlbAppr(k)=det(BysCrlbAppr(:,:,k));
end


figure
semilogy(SNRdB,DetMlCov); hold on;
semilogy(SNRdB,DetMlCrlbNum)
semilogy(SNRdB,DetBysCov);
semilogy(SNRdB,DetBysCrlbNum);
xlabel 'SNR (dB)'; ylabel 'CRLB and Error Cov. Det. (mV^2.Sec^4)';
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
title 'Det. of Error Cov. and CRLB matrices for Gaus. parameters'
saveas(gcf,'DetOneGauss_BadInitialized.fig')
saveas(gcf,'DetOneGauss_BadInitialized.eps','epsc')
% saveas(gcf,'DetOneGauss_GoodInitialized.fig')
% saveas(gcf,'DetOneGauss_GoodInitialized.eps','epsc')


figure
semilogy(SNRdB,DetMlCrlbNum); hold on;
semilogy(SNRdB,DetMlCrlbAppr)
semilogy(SNRdB,DetBysCrlbNum);
semilogy(SNRdB,DetBysCrlbAppr);
xlabel 'SNR (dB)'; ylabel 'CRLB Det. (mV^2.Sec^4)';
legend('ML Num. CRLB','ML Appr. CRLB',' BYS Num. CRLB','BYS Appr. CRLB')
title 'Det. of Numeric CRLB vs Approximated CRLB matrices for a Gaus. parameters'
saveas(gcf,'NumApprDetCrlbOneGauss.fig')
saveas(gcf,'NumApprDetCrlbOneGauss.eps','epsc')







function g=gauss(t,pp)

    g=pp(1)*exp(-((t-pp(3)).^2)./(2*pp(2)^2));
end
