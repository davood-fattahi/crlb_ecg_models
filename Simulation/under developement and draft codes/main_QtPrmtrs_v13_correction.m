clear 
close all
clc

Addrss='Backup and Results';
Files=dir([Addrss '\*TempQt*.mat']); 
NCases=length(Files);

for ii=1:NCases
    load([Addrss '\' Files(ii).name]); % address of the dataset   
   
    %% Stage 2- estimating the parameters from noisy signals
    h = waitbar(0,['Case No. ' num2str(ii) ', Estimating parameters for noisy signals, please wait ...']);

    
    %%% pre allocating ...
    BysParamsQ=zeros(NumGaus*3, NumBeats, NumRuns, length(SNR));
    BysErQ=zeros(size(BysParamsQ));     
    BysParamsT=zeros(NumGaus*3, NumBeats, NumRuns, length(SNR));
    BysErT=zeros(size(BysParamsT));     
    BysCovQ=zeros(3,3,length(SNR));
    DetBysCovQ=zeros(size(SNR)); 
    BysCovT=zeros(3,3,length(SNR));
    DetBysCovT=zeros(size(SNR)); 
 
 
   
    %%% 
    for k=1:length(SNR) % for each SNR
        rng(k+i); % fix the seed for noise generating
        signalNQ=signal+randn(length(signal),1).*sqrt(VarNoiseAddQ(k)); % adding noise to the signal
        rng(k+i+1); % fix the seed for noise generating
        signalNT=signal+randn(length(signal),1).*sqrt(VarNoiseAddT(k)); % adding noise to the signal
        for m=1:NumRuns % for each run
            waitbar(((k-1)*NumRuns+m)/(length(SNR)*NumRuns)) % wait bar ...
            jj=0; % initializing
            for j=JJ
                jj=jj+1; % number of beats counter
                SegNQ=signalNQ(QSegOn(j):QSegOff(j)); % SOI
                tt=tm(QSegOn(j):QSegOff(j))-Rt(j,2); % SOI time stamp
                p0=[signal(Q(j,2)) (Qt(j,3)-Qt(j,1))/5 Qt(j,2)]; % initial points
                 %%% for Bayesian, fit gaussians on each noisy Q-segs
                BysParamsQ(:,jj,m,k)=GausFit(tt,SegNQ,p0,eval(Qplb),eval(Qpub),PrMeanQ,PrCovQ,VarSigQ/SNR(k),optionsBys);
                
                SegNT=signalNT(TSegOn(j):TSegOff(j)); % SOI
                tt=tm(TSegOn(j):TSegOff(j))-Rt(j,2); % SOI time stamp
                p0=[signal(T(j,2)) (Tt(j,3)-Tt(j,1))/5 Tt(j,2)]; % initial points
                %%% for Bayesian, fit gaussians on each noisy T-segs
                BysParamsT(:,jj,m,k)=GausFit(tt,SegNT,p0,eval(Tplb),eval(Tpub),PrMeanT,PrCovT,VarSigT/SNR(k),optionsBys);
            end
            BysErQ(:,:,m,k)=BysParamsQ(:,:,m,k)-PrParamsQ(JJ,:)'; % Bys Error from prior
            BysErT(:,:,m,k)=BysParamsT(:,:,m,k)-PrParamsT(JJ,:)'; % Bys Error from prior         
        end
        BysCovQ(:,:,k)=cov(reshape(BysErQ(:,:,:,k),3,[])'); % Bys Error covariance matrix
        DetBysCovQ(k)=det(BysCovQ(:,:,k)); % determinant of Bys Error covariance matrix
        
        BysCovT(:,:,k)=cov(reshape(BysErT(:,:,:,k),3,[])'); % Bys Error covariance matrix
        DetBysCovT(k)=det(BysCovT(:,:,k)); % determinant of Bys Error covariance matrix
    end
    close(h)


%     save(['Backup and Results\TempQtPrmtrs_' files(i).name(1:end-4)  '.mat'], 'PrParamsQ', 'PrParamsT', 'MlParamsQ', 'BysParamsQ', 'MlParamsT', 'BysParamsT', 'VarNoiseAddQ', 'VarNoiseAddT', 'VarSigQ', 'VarSigT', 'SNR', 'MlCovQ', 'DetMlCovQ', 'BysCovQ', 'DetBysCovQ', 'MlCovT', 'DetMlCovT', 'BysCovT', 'DetBysCovT', 'MlMseQ', 'BysMseQ', 'MlMseT', 'BysMseT', 'MlCrlbQ', 'BysCrlbQ', 'MlCrlbT', 'BysCrlbT', 'MlCrlbQ_AvrgdFI', 'MlCrlbT_AvrgdFI',  'MlFIQ', 'MlFIT', 'DetMlCrlbQ', 'DetBysCrlbQ', 'DetMlCrlbT', 'DetBysCrlbT')    
    save(['Backup and Results\TempQtPrmtrs_' files(i).name(1:end-4)  '_Crct.mat']);
end
% 
% figure
% semilogy(SNRdB,DetMlCrlbNumQ)
% hold on
% semilogy(SNRdB,DetMlCrlbApprQ)
% xlabel 'SNR (dB)'; ylabel 'CRLB Det. (mV^2.Sec^4)';
% legend('Num. CRLB','Appr. CRLB')
% title 'Det. of Numeric CRLB vs. Approximated CRLB for Q wave parameters'
% saveas(gcf,'NumApprCrlbDetQ.fig')
% saveas(gcf,'NumApprCrlbDetQ.eps','epsc')
% 
% 
% figure
% semilogy(SNRdB,DetMlCrlbNumT)
% hold on
% semilogy(SNRdB,DetMlCrlbApprT)
% xlabel 'SNR (dB)'; ylabel 'CRLB Det. (mV^2.Sec^4)';
% legend('Num. CRLB','Appr. CRLB')
% title 'Det. of Numeric CRLB vs. Approximated CRLB for T wave parameters'
% saveas(gcf,'NumApprCrlbDetT.fig')
% saveas(gcf,'NumApprCrlbDetT.eps','epsc')
% 
% 
% 
% figure
% semilogy(SNRdB,DetMlCovQ); hold on;
% semilogy(SNRdB,DetMlCrlbNumQ)
% semilogy(SNRdB,DetBysCovQ);
% semilogy(SNRdB,DetBysCrlbNumQ);
% xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% % title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
% saveas(gcf,'DetQ.fig')
% saveas(gcf,'DetQ.eps','epsc')
% 
% 
% figure
% semilogy(SNRdB,DetMlCovT); hold on;
% semilogy(SNRdB,DetMlCrlbNumT)
% semilogy(SNRdB,DetBysCovT);
% semilogy(SNRdB,DetBysCrlbNumT);
% xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% % title 'Determinants of Error Cov. and CRLB matrices for T wave Gaus. parameters'
% saveas(gcf,'DetT.fig')
% saveas(gcf,'DetT.eps','epsc')
