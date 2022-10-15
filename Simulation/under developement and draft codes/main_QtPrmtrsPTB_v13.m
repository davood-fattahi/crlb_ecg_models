clc
clear
close all


%% the constants
N=[];
N0=1;
Anttr='q1c'; 
Chnl = 2; % 2->ii, 8->v2, 11->v5
HrTr= 120; % bpm, threshold for heart rate
L1=60; % ms, short st length
L2=80; % ms, long st length
SNRdB=(-20:5:40); % vector
SNR=10.^(SNRdB/10);
NumGaus=1; 
w1=.75;
w2=.9;
NumBeats=100;
NumRuns=1;
AllBeats=1000;


% optionsMl = struct('Display','final-detailed' ,'MaxFunctionEvaluations',10000, 'MaxIterations', 10000, 'OptimalityTolerance', 1e-20, 'StepTolerance', 1e-20);
% optionsBys = struct('Display','final-detailed' ,'MaxFunctionEvaluations',10000, 'MaxIterations', 10000, 'ConstraintTolerance', 1e-20, 'OptimalityTolerance', 1e-20, 'StepTolerance', 1e-20);
optionsMl=struct;
optionsBys=struct;

Qsegon=[1 0]; % determining the SOI begining time (ms); Qsegon=[alpha beta] -> Qsegon = time(alpha) + beta; alpha =1,2,3 respectively for begining, center, and end of the wave.
Qsegoff=[3 0]; % determining the SOI end time (ms); Qsegoff=[alpha beta] -> Qsegoff = time(alpha) + beta; alpha =1,2,3 respectively for begining, center, and end of the wave.
Qplb='[-inf 0 Qt(j,1)]';
Qpub='[0 Qt(j,3)-Qt(j,1) Qt(j,3)]';

Tsegon=[1 0]; % determining the SOI begining time (ms); Tsegon=[alpha beta] -> Tsegon = time(alpha) + beta; alpha =1,2,3 respectively for begining, center, and end of the wave.
Tsegoff=[3 0]; % determining the SOI end time (ms); Tsegoff=[alpha beta] -> Tsegoff = time(alpha) + beta; alpha =1,2,3 respectively for begining, center, and end of the wave.
Tplb='[0 0 Tt(j,1)]';
Tpub='[inf Tt(j,3)-Tt(j,1) Tt(j,3)]';

% for example, in order to select a segment from 20 ms befor the center of T wave to the end of it:
% Tsegon=[2 -20]; 
% Tsegoff=[3 0];


%% load the case names in the directory
NC=importdata('..\mcode\database\ptb-diagnostic-ecg-database-1.0.0\CONTROLS'); % load the normal cases
NumCases=length(NC);
DA='database\ptb-diagnostic-ecg-database-1.0.0';
for i=1:NumCases
    %% load the record
    oldFolder=cd('..\mcode');
    AdrsNm=[DA  '\' NC{i}];
    [signal,Fs,tm]=rdsamp(AdrsNm,Chnl,N,N0); % 2->ii, 8->v2, 11->v5
    [annatr,anntypeatr,subtypeatr,chanatr,numatr,commentsatr]=rdann(AdrsNm, 'hea', [], N, N0);

    cd(oldFolder)


    % extract the RR-ints
    Rp = PeakDetection(signal,1/Fs); Rp=find(Rp);


    if isempty(AllBeats) || (AllBeats >length(Rp)-2)
        NumAllBeats=length(Rp)-2;
    else
        NumAllBeats=AllBeats;
    end
    
    %% baseline wandering removal,
    signal=signal-(BaseLine1(BaseLine1(signal', round(w1*Fs), 'md'), round(w2*Fs), 'mn'))';
    
    % extract approx. fiducial points
    [P, Q, R, S, T] =ecgfidapprox2(signal, Rp);
    P([1 end],:)=[]; Q([1 end],:)=[]; R([1 end],:)=[]; S([1 end],:)=[]; T([1 end],:)=[];
    Rt=tm(R); Qt=tm(Q)-Rt(:,2); Tt=tm(T)-Rt(:,2);
    QSegOn=Q(:,Qsegon(1))+floor(Qsegon(2)*Fs/1000);
    QSegOff=Q(:,Qsegoff(1))+floor(Qsegoff(2)*Fs/1000);
    TSegOn=T(:,Tsegon(1))+floor(Tsegon(2)*Fs/1000);
    TSegOff=T(:,Tsegoff(1))+floor(Tsegoff(2)*Fs/1000);
    
    PlotFP(tm,signal, P(:), Q(:), R(:), S(:), T(:))
    saveas(gcf,'ApprFP.fig')


    %% computing the prior
    % Predefine ...
    PrParamsQ=zeros(NumAllBeats,NumGaus*3);
    PrParamsT=zeros(NumAllBeats,NumGaus*3);

    
    % wait bar
    h = waitbar(0,['Case No. ' num2str(i) ', Estimating parameters for clean signals, please wait ...']);
    LQ=0; LT=0; SegEvsQ=[]; SegsEvT=[]; 
    SegsQ=nan(sum((QSegOff(1:NumAllBeats)-QSegOn(1:NumAllBeats)+1)),1);
    SegsT=nan(sum((TSegOff(1:NumAllBeats)-TSegOn(1:NumAllBeats)+1)),1);

    for j=1:NumAllBeats
        waitbar(j/NumAllBeats)
        %% Q wave model fitting
        SegQ=signal(QSegOn(j):QSegOff(j));
        SegsQ(LQ+1:LQ+length(SegQ))=SegQ(:);
        LQ=LQ+length(SegQ);
        tt=tm(QSegOn(j):QSegOff(j))-Rt(j,2);
        
        % fit gssns on each clean Q-wave
        p0=[signal(Q(j,2)) (Qt(j,3)-Qt(j,1))/5 Qt(j,2)];
        PrParamsQ(j,:)=GausFit(tt,SegQ,p0,eval(Qplb),eval(Qpub),optionsMl);
        
        % polt
        SegEvQ=GausVal(tt,PrParamsQ(j,:));
        SegEvsQ=[SegEvsQ; SegEvQ];        
        hold on;
        plot(tt+Rt(j,2),SegEvQ,'r-')
       
        %% T wave model fitting      
        SegT=signal(TSegOn(j):TSegOff(j));
        SegsT(LT+1:LT+length(SegT))=SegT(:);
        LT=LT+length(SegT);
        tt=tm(TSegOn(j):TSegOff(j))-Rt(j,2);
        
        % fit gssns on each clean T-halfwave
        p0=[signal(T(j,2)) (Tt(j,3)-Tt(j,1))/5 Tt(j,2)];
        PrParamsT(j,:)=GausFit(tt,SegT,p0,eval(Tplb),eval(Tpub),optionsMl);
        
        % plot
        tt=tm(T(j,1):T(j,3))-Rt(j,2);
        SegEvT=GausVal(tt(:),PrParamsT(j,:));
        SegsEvT=[SegsEvT; SegEvT];     
        hold on
        plot(tt(:)+Rt(j,2),SegEvT,'r-')

    end
    close(h)

%     figure
%     histogram(QPrParams(:,1),100);
%     figure
%     histogram(QPrParams(:,2));
%     figure
%     histogram(QPrParams(:,3));
%     
%     figure
%     histogram(TPrParams(:,1),100);
%     figure
%     histogram(TPrParams(:,2));
%     figure
%     histogram(TPrParams(:,3));
%     
%     %% Error histogram
%     QMdlErr=QSegs-QSegEvs;
%     TMdlErr=TSegs-TSegEvs;
% 
%     figure; histogram(QMdlErr(:)); 
%     figure; histogram(TMdlErr(:)); legend('Model Error of Q', 'Model Error of T' ); 
%     
    %%% estimate the prior for params (cov matrix of params)
    PrMeanQ=mean(PrParamsQ,1);
    PrCovQ=cov(PrParamsQ);
    VarSigQ=var(SegsQ(:)); % variance of the Q signal
    PrMeanT=mean(PrParamsT,1);
    PrCovT=cov(PrParamsT);
    VarSigT=var(SegsT(:)); % variance of the T signal
    
    %% add noise to the clean Q-seg
    h = waitbar(0,['Case No. ' num2str(i) ', Estimating parameters for noisy signals, please wait ...']);
    if isempty(NumBeats)
        NumBeats=NumAllBeats;
    end
    if NumBeats > NumAllBeats
        NumBeats=NumAllBeats;
    end
    MlParamsQ=zeros(NumGaus*3, NumBeats, NumRuns, length(SNR));
    MlErQ=zeros(size(MlParamsQ));     
    
    BysParamsQ=zeros(NumGaus*3, NumBeats, NumRuns, length(SNR));
    BysErQ=zeros(size(BysParamsQ));     
    
    MlParamsT=zeros(NumGaus*3, NumBeats, NumRuns, length(SNR)); 
    MlErT=zeros(size(MlParamsT)); 
    
    BysParamsT=zeros(NumGaus*3, NumBeats, NumRuns, length(SNR));
    BysErT=zeros(size(BysParamsT));     
   
    VarNoiseAddQ=zeros(size(SNR)); 
    VarNoiseAddT=zeros(size(SNR)); 
    
    MlCovQ=zeros(3,3,length(SNR));
    BysCovQ=zeros(3,3,length(SNR));
   
    DetMlCovQ=zeros(size(SNR)); 
    DetBysCovQ=zeros(size(SNR)); 
   
    MlCovT=zeros(3,3,length(SNR));
    BysCovT=zeros(3,3,length(SNR));

    DetMlCovT=zeros(size(SNR)); 
    DetBysCovT=zeros(size(SNR)); 
 
 
        
    rng(i)
    JJ=randperm(NumAllBeats,NumBeats);
    
    for k=1:length(SNR)
        VarNoiseAddQ(k)=VarSigQ/SNR(k); 
        VarNoiseAddT(k)=VarSigT/SNR(k); 
        rng(k+i); signalNQ=signal+randn(length(signal),1).*sqrt(VarNoiseAddQ(k));
        rng(k+i+1); signalNT=signal+randn(length(signal),1).*sqrt(VarNoiseAddT(k));
        for m=1:NumRuns
            waitbar(((k-1)*NumRuns+m)/(length(SNR)*NumRuns))
            jj=0;
            for j=JJ
                jj=jj+1;
                SegNQ=signalNQ(QSegOn(j):QSegOff(j));
                tt=tm(QSegOn(j):QSegOff(j))-Rt(j,2);
                p0=[signal(Q(j,2)) (Qt(j,3)-Qt(j,1))/5 Qt(j,2)];
                %%% for ML, fit gaussians on each noisy Q-segs
                MlParamsQ(:,jj,m,k)=GausFit(tt,SegNQ,p0,eval(Qplb),eval(Qpub),optionsMl);
                %%% for Bayesian, fit gaussians on each noisy Q-segs
                BysParamsQ(:,jj,m,k)=GausFit(tt,SegNQ,p0,eval(Qplb),eval(Qpub),PrMeanQ,PrCovQ,VarSigQ/SNR(k),optionsBys);
                
                
                SegNT=signalNT(TSegOn(j):TSegOff(j));
                tt=tm(TSegOn(j):TSegOff(j))-Rt(j,2);
                p0=[signal(T(j,2)) (Tt(j,3)-Tt(j,1))/5 Tt(j,2)];
                %%% for ML, fit gaussians on each noisy T-segs
                MlParamsT(:,jj,m,k)=GausFit(tt,SegNT,p0,eval(Tplb),eval(Tpub),optionsMl);
                %%% for Bayesian, fit gaussians on each noisy T-segs
                BysParamsT(:,jj,m,k)=GausFit(tt,SegNT,p0,eval(Tplb),eval(Tpub),PrMeanT,PrCovT,VarSigT/SNR(k),optionsBys);
            end
            MlErQ(:,:,m,k)=MlParamsQ(:,:,m,k)-PrParamsQ(JJ,:)';
            BysErQ(:,:,m,k)=BysParamsQ(:,:,m,k)-PrParamsQ(JJ,:)';
            MlErT(:,:,m,k)=MlParamsT(:,:,m,k)-PrParamsT(JJ,:)';
            BysErT(:,:,m,k)=BysParamsT(:,:,m,k)-PrParamsT(JJ,:)';            
        end
        MlCovQ(:,:,k)=cov(reshape(MlErQ(:,:,:,k),3,[])');
        DetMlCovQ(k)=det(MlCovQ(:,:,k));
        BysCovQ(:,:,k)=cov(reshape(BysErQ(:,:,:,k),3,[])');
        DetBysCovQ(k)=det(BysCovQ(:,:,k));
        
        MlCovT(:,:,k)=cov(reshape(MlErT(:,:,:,k),3,[])');
        DetMlCovT(k)=det(MlCovT(:,:,k));
        BysCovT(:,:,k)=cov(reshape(BysErT(:,:,:,k),3,[])');
        DetBysCovT(k)=det(BysCovT(:,:,k));
    end
    close(h)
    MlMseQ=squeeze(mean((MlParamsQ-PrParamsQ(JJ,:)').^2, [2 3]));
    BysMseQ=squeeze(mean((BysParamsQ-PrParamsQ(JJ,:)').^2, [2 3]));
    MlMseT=squeeze(mean((MlParamsT-PrParamsT(JJ,:)').^2, [2 3]));
    BysMseT=squeeze(mean((BysParamsT-PrParamsT(JJ,:)').^2, [2 3]));

%% CRLB calculation
    MlCrlbApprQ=zeros(3, 3, NumAllBeats, length(SNRdB));
    MlFIApprQ=zeros(3, 3, NumAllBeats, length(SNRdB));
    BysCrlbApprQ=zeros(3, 3, length(SNRdB));
    MlCrlbApprT=zeros(3, 3, NumAllBeats, length(SNRdB));
    MlFIApprT=zeros(3, 3, NumAllBeats, length(SNRdB));
    BysCrlbApprT=zeros(3, 3, length(SNRdB));
    MlCrlbApprQ_AvrgdFI=zeros(3, 3, length(SNRdB));
    MlCrlbApprT_AvrgdFI=zeros(3, 3, length(SNRdB));
    DetMlCrlbApprQ=zeros(length(SNRdB),1);
    DetBysCrlbApprQ=zeros(length(SNRdB),1);
    DetMlCrlbApprT=zeros(length(SNRdB),1);
    DetBysCrlbApprT=zeros(length(SNRdB),1);
    
    MlCrlbNumQ=zeros(3, 3, NumAllBeats, length(SNRdB));
    MlFINumQ=zeros(3, 3, NumAllBeats, length(SNRdB));
    BysCrlbNumQ=zeros(3, 3, length(SNRdB));
    MlCrlbNumT=zeros(3, 3, NumAllBeats, length(SNRdB));
    MlFINumT=zeros(3, 3, NumAllBeats, length(SNRdB));
    BysCrlbNumT=zeros(3, 3, length(SNRdB));
    MlCrlbNumQ_AvrgdFI=zeros(3, 3, length(SNRdB));
    MlCrlbNumT_AvrgdFI=zeros(3, 3, length(SNRdB));
    DetMlCrlbNumQ=zeros(length(SNRdB),1);
    DetBysCrlbNumQ=zeros(length(SNRdB),1);
    DetMlCrlbNumT=zeros(length(SNRdB),1);
    DetBysCrlbNumT=zeros(length(SNRdB),1);
    
    for k=1:length(SNRdB)
        for j=1:NumAllBeats
            [MlCrlbNumQ(:,:,j,k), MlFINumQ(:,:,j,k)]=GaussCrlbNumeric(tm(QSegOn(j):QSegOff(j))-Rt(j,2), PrParamsQ(j,:), VarNoiseAddQ(k), 1 );
            [MlCrlbNumT(:,:,j,k), MlFINumT(:,:,j,k)]=GaussCrlbNumeric(tm(TSegOn(j):TSegOff(j))-Rt(j,2), PrParamsT(j,:), VarNoiseAddT(k), 1 );
            
            [MlCrlbApprQ(:,:,j,k), MlFIApprQ(:,:,j,k)]=GaussCRLB(PrParamsQ(j,1), PrParamsQ(j,2), Fs, VarNoiseAddQ(k) );
            [MlCrlbApprT(:,:,j,k), MlFIApprT(:,:,j,k)]=GaussCRLB(PrParamsT(j,1), PrParamsT(j,2), Fs, VarNoiseAddT(k) );
%             [crlb, fi]=GaussCRLB(PrParamsT(j,1), PrParamsT(j,2), Fs, VarNoiseAddT(k) );
%             MlCrlbApprT(:,:,j,k)=2*crlb;  MlFIApprT(:,:,j,k)=.5*fi; % since we consider half of T wave        
        end
        MlCrlbNumQ_AvrgdFI(:,:,k)=inv(mean(MlFINumQ(:,:,:,k),3));
        BysCrlbNumQ(:,:,k)=inv(mean(MlFINumQ(:,:,:,k),3)+inv(PrCovQ));

        MlCrlbNumT_AvrgdFI(:,:,k)=inv(mean(MlFINumT(:,:,:,k),3));
        BysCrlbNumT(:,:,k)=inv(mean(MlFINumT(:,:,:,k),3)+inv(PrCovT));

        DetMlCrlbNumQ(k)=det(MlCrlbNumQ_AvrgdFI(:,:,k));
        DetBysCrlbNumQ(k)=det(BysCrlbNumQ(:,:,k));
%         EigBysCrlbQ(:,k)=eig(BysCovQ(:,:,k)-BysCrlbQ(:,:,k));
%         EigMlCrlbQ(:,k)=eig(MlCovQ(:,:,k)-MlCrlbQ_AvrgdFI(:,:,k));   

        DetMlCrlbNumT(k)=det(MlCrlbNumT_AvrgdFI(:,:,k));
        DetBysCrlbNumT(k)=det(BysCrlbNumT(:,:,k));
%         EigBysCrlbT(:,k)=eig(BysCovT(:,:,k)-BysCrlbT(:,:,k));
%         EigMlCrlbQ(:,k)=eig(MlCovT(:,:,k)-MlCrlbT_AvrgdFI(:,:,k));

        MlCrlbApprQ_AvrgdFI(:,:,k)=inv(mean(MlFIApprQ(:,:,:,k),3));
        BysCrlbApprQ(:,:,k)=inv(mean(MlFIApprQ(:,:,:,k),3)+inv(PrCovQ));

        MlCrlbApprT_AvrgdFI(:,:,k)=inv(mean(MlFIApprT(:,:,:,k),3));
        BysCrlbApprT(:,:,k)=inv(mean(MlFIApprT(:,:,:,k),3)+inv(PrCovT));

        DetMlCrlbApprQ(k)=det(MlCrlbApprQ_AvrgdFI(:,:,k));
        DetBysCrlbApprQ(k)=det(BysCrlbApprQ(:,:,k));
%         EigBysCrlbQ(:,k)=eig(BysCovQ(:,:,k)-BysCrlbQ(:,:,k));
%         EigMlCrlbQ(:,k)=eig(MlCovQ(:,:,k)-MlCrlbQ_AvrgdFI(:,:,k));   

        DetMlCrlbApprT(k)=det(MlCrlbApprT_AvrgdFI(:,:,k));
        DetBysCrlbApprT(k)=det(BysCrlbApprT(:,:,k));
%         EigBysCrlbT(:,k)=eig(BysCovT(:,:,k)-BysCrlbT(:,:,k));
%         EigMlCrlbQ(:,k)=eig(MlCovT(:,:,k)-MlCrlbT_AvrgdFI(:,:,k));
    end

%     save(['Backup and Results\TempQtPrmtrs_' files(i).name(1:end-4)  '.mat'], 'PrParamsQ', 'PrParamsT', 'MlParamsQ', 'BysParamsQ', 'MlParamsT', 'BysParamsT', 'VarNoiseAddQ', 'VarNoiseAddT', 'VarSigQ', 'VarSigT', 'SNR', 'MlCovQ', 'DetMlCovQ', 'BysCovQ', 'DetBysCovQ', 'MlCovT', 'DetMlCovT', 'BysCovT', 'DetBysCovT', 'MlMseQ', 'BysMseQ', 'MlMseT', 'BysMseT', 'MlCrlbQ', 'BysCrlbQ', 'MlCrlbT', 'BysCrlbT', 'MlCrlbQ_AvrgdFI', 'MlCrlbT_AvrgdFI',  'MlFIQ', 'MlFIT', 'DetMlCrlbQ', 'DetBysCrlbQ', 'DetMlCrlbT', 'DetBysCrlbT')    
    save(['Backup and Results\TempQt2Prmtrs_' NC{i}(end-7:end)  '.mat']);
end

figure
semilogy(SNRdB,DetMlCrlbNumQ)
hold on
semilogy(SNRdB,DetMlCrlbApprQ)
xlabel 'SNR (dB)'; ylabel 'CRLB Det. (mV^2.Sec^4)';
legend('Num. CRLB','Appr. CRLB')
title 'Det. of Numeric CRLB vs. Approximated CRLB for Q wave parameters'
saveas(gcf,'NumApprCrlbDetQ.fig')
saveas(gcf,'NumApprCrlbDetQ.eps','epsc')


figure
semilogy(SNRdB,DetMlCrlbNumT)
hold on
semilogy(SNRdB,DetMlCrlbApprT)
xlabel 'SNR (dB)'; ylabel 'CRLB Det. (mV^2.Sec^4)';
legend('Num. CRLB','Appr. CRLB')
title 'Det. of Numeric CRLB vs. Approximated CRLB for T wave parameters'
saveas(gcf,'NumApprCrlbDetT.fig')
saveas(gcf,'NumApprCrlbDetT.eps','epsc')



figure
semilogy(SNRdB,DetMlCovQ); hold on;
semilogy(SNRdB,DetMlCrlbNumQ)
semilogy(SNRdB,DetBysCovQ);
semilogy(SNRdB,DetBysCrlbNumQ);
xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'DetQ.fig')
saveas(gcf,'DetQ.eps','epsc')


figure
semilogy(SNRdB,DetMlCovT); hold on;
semilogy(SNRdB,DetMlCrlbNumT)
semilogy(SNRdB,DetBysCovT);
semilogy(SNRdB,DetBysCrlbNumT);
xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'Determinants of Error Cov. and CRLB matrices for T wave Gaus. parameters'
saveas(gcf,'DetT.fig')
saveas(gcf,'DetT.eps','epsc')
