clc
clear
close all


%% the constants
N=[]; % the end sample of the signal loaded by rdsamp
N0=1; % the start sample of the signal loaded by rdsamp
Anttr='q1c'; % the annoattor name used in R peak detection
Chnl = 1; % ECG channel 
HrTr= 70; % bpm, threshold for heart rate, used in r peak detector
SNRdB=(-20:5:40); % vector of SNRs in decibel
SNR=10.^(SNRdB/10); % SNRs
NumGaus=1; % number of Gaussians utilized in modeling ecg segments
w1=.75; % window length of the first median filter used in base line removing
w2=.9; % window length of the second median filter used in base line removing
NumBeats=1000; % number of the beats envolved in noisy signals' parameter estimation stage
NumRuns=3; % number of repeats in noisy signals' parameter estimation stage
AllBeats=[]; % number of the beats envolved in prior parameter estimation stage (clean signals)

%%% option structure for non linear least squre optimization, respectively
%%% in ML estimation and Bayesian estimation:
% optionsMl = struct('Display','final-detailed' ,'MaxFunctionEvaluations',10000, 'MaxIterations', 10000, 'OptimalityTolerance', 1e-20, 'StepTolerance', 1e-20);
% optionsBys = struct('Display','final-detailed' ,'MaxFunctionEvaluations',10000, 'MaxIterations', 10000, 'ConstraintTolerance', 1e-20, 'OptimalityTolerance', 1e-20, 'StepTolerance', 1e-20);
optionsMl=struct;
optionsBys=struct;

Qsegon=[1 0]; % determining the SOI begining time (ms); Qsegon=[alpha beta] -> Qsegon = time(alpha) + beta; alpha =1,2,3 respectively for begining, center, and end of the wave.
Qsegoff=[3 0]; % determining the SOI end time (ms); Qsegoff=[alpha beta] -> Qsegoff = time(alpha) + beta; alpha =1,2,3 respectively for begining, center, and end of the wave.
Qplb='[-inf 0 Qt(j,1)]'; % optimization conditions: lower bound for the Q wave parameters
Qpub='[inf Qt(j,3)-Qt(j,1) Qt(j,3)]'; % optimization conditions: upper bound for the Q wave parameters

Tsegon=[2 -20]; % determining the SOI begining time (ms); Tsegon=[alpha beta] -> Tsegon = time(alpha) + beta; alpha =1,2,3 respectively for begining, center, and end of the wave.
Tsegoff=[3 0]; % determining the SOI end time (ms); Tsegoff=[alpha beta] -> Tsegoff = time(alpha) + beta; alpha =1,2,3 respectively for begining, center, and end of the wave.
Tplb='[-inf 0 Tt(j,1)]';% optimization conditions: lower bound for the T wave parameters
Tpub='[inf Tt(j,3)-Tt(j,1) Tt(j,3)]';% optimization conditions: upper bound for the T wave parameters


% for example, in order to select a segment from 20 ms befor the center of T wave to the end of it:
% Tsegon=[2 -20]; 
% Tsegoff=[3 0];


%% load the case names in the directory
DA='..\mcode\database\qt-database-1.0.0';
files=dir([DA '\*.dat']); 
NumCases=length(files);
DA='database\qt-database-1.0.0';

for i=1:NumCases
    %% load the record
    oldFolder=cd('..\mcode'); % go to the wfdb toolbox root
    AdrsNm=[DA  '\' files(i).name(1:end-4)]; % address of the dataset
    [signal,Fs,tm]=rdsamp(AdrsNm,[],N,N0); % loading the data
    
    %% extract the R peaks using annotation
%     try
%     [RR,Rp]=ann2rr(AdrsNm,'atr',N, N0);
%     catch ME
%         warning(['Case No. ' num2str(i) ', ' files(i).name(1:end-4) ' is skiped. See the error bellow:' newline ME.message]);
%         cd(oldFolder)
%         continue
%     end
%     
 

    cd(oldFolder) % go back to the main folder
    
    %% baseline wandering removal,
    signal=signal-(BaseLine1(BaseLine1(signal', round(w1*Fs), 'md'), round(w2*Fs), 'mn'))';
 

    %% extract the R peaks using peak detector
    Rp = PeakDetection(mean(signal,2),1/Fs); Rp=find(Rp);
  
    %% adjusting the AllBeats by the detected R peaks:
    if isempty(AllBeats) || (AllBeats >length(Rp)-2)  %#ok<BDSCI>
        NumAllBeats=length(Rp)-2;
    else
        NumAllBeats=AllBeats;
    end
    
   
    %% extract approx. fiducial points
    signal=signal(:,Chnl);
    
    [~, Q, R, ~, T] =ecgfidapprox2(signal, Rp);
    Q([1 end],:)=[]; R([1 end],:)=[]; T([1 end],:)=[];
    Rt=tm(R); Qt=tm(Q)-Rt(:,2); Tt=tm(T)-Rt(:,2); % getting time stamps of the fid points
    
    %%% adjusting the start and the end of SOI
    QSegOn=Q(:,Qsegon(1))+floor(Qsegon(2)*Fs/1000); 
    QSegOff=Q(:,Qsegoff(1))+floor(Qsegoff(2)*Fs/1000); 
    TSegOn=T(:,Tsegon(1))+floor(Tsegon(2)*Fs/1000);
    TSegOff=T(:,Tsegoff(1))+floor(Tsegoff(2)*Fs/1000);
    
    %% plot the fid points on the signal
%     figure
%     plot(tm,signal); hold on
%     plot(tm(Q(~isnan(Q))),signal(Q(~isnan(Q))),'*')
%     plot(tm(Rp(~isnan(Rp))),signal(Rp(~isnan(Rp))),'*')
%     plot(tm(T(~isnan(T))),signal(T(~isnan(T))),'*')    
%     saveas(gcf,'ApprFP.fig')


    %% stage 1- computing the prior
    % Pre-allocating ...
    PrParamsQ=zeros(NumAllBeats,NumGaus*3);
    PrParamsT=zeros(NumAllBeats,NumGaus*3);
    LQ=0; LT=0; SegEvsQ=[]; SegsEvT=[]; 
    SegsQ=nan(sum((QSegOff(1:NumAllBeats)-QSegOn(1:NumAllBeats)+1)),1);
    SegsT=nan(sum((TSegOff(1:NumAllBeats)-TSegOn(1:NumAllBeats)+1)),1);

    
    % wait bar
    h = waitbar(0,['Case No. ' num2str(i) ', Estimating parameters for clean signals, please wait ...']);

    for j=1:NumAllBeats  % for each beat ...
        waitbar(j/NumAllBeats)
        %% Q wave model fitting
        SegQ=signal(QSegOn(j):QSegOff(j)); % Segment Of Intrest (SOI)
        SegsQ(LQ+1:LQ+length(SegQ))=SegQ(:); % All the SOIs
        LQ=LQ+length(SegQ); % All the SOIs length
        tt=tm(QSegOn(j):QSegOff(j))-Rt(j,2); % SOI time stamp
        
        % fit gssns on each clean Q-wave
        p0=[signal(Q(j,2)) (Qt(j,3)-Qt(j,1))/5 Qt(j,2)]; % initial condition
        PrParamsQ(j,:)=GausFit(tt,SegQ,p0,eval(Qplb),eval(Qpub),optionsMl);
        
%         % polt the evaluated gaussians on the signal
%         SegEvQ=GausVal(tt,PrParamsQ(j,:));
%         hold on;
%         plot(tt+Rt(j,2),SegEvQ,'r-')
%        
        %% T wave model fitting      
        SegT=signal(TSegOn(j):TSegOff(j)); % Segment Of Intrest (SOI)
        SegsT(LT+1:LT+length(SegT))=SegT(:); % All the SOIs
        LT=LT+length(SegT); % All the SOIs length
        tt=tm(TSegOn(j):TSegOff(j))-Rt(j,2); % SOI time stamp
        
        % fit gssns on each clean T-halfwave
        p0=[signal(T(j,2)) (Tt(j,3)-Tt(j,1))/5 Tt(j,2)]; % initial condition
        PrParamsT(j,:)=GausFit(tt,SegT,p0,eval(Tplb),eval(Tpub),optionsMl);
        
%         % polt the evaluated gaussians on the signal
%         tt=tm(T(j,1):T(j,3))-Rt(j,2);
%         SegEvT=GausVal(tt(:),PrParamsT(j,:));
%         hold on
%         plot(tt(:)+Rt(j,2),SegEvT,'r-')

    end
    close(h)
    
%     %% plot histogram of the prior parameters ...
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
    %% estimate the prior for the params
    PrMeanQ=mean(PrParamsQ,1);
    PrCovQ=cov(PrParamsQ);
    VarSigQ=var(SegsQ(:)); % variance of the Q signal
    PrMeanT=mean(PrParamsT,1);
    PrCovT=cov(PrParamsT);
    VarSigT=var(SegsT(:)); % variance of the T signal
    
    %% Stage 2- estimating the parameters from noisy signals
    h = waitbar(0,['Case No. ' num2str(i) ', Estimating parameters for noisy signals, please wait ...']);
    
    %%% adjusting the number of noisy beats
    if isempty(NumBeats)
        NumBeats=NumAllBeats;
    end
    if NumBeats > NumAllBeats
        NumBeats=NumAllBeats;
    end
    
    %%% pre allocating ...
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
 
 
    %%% select the noisy beats randomly
    rng(i) % fixing the seed for random permuting
    JJ=randperm(NumAllBeats,NumBeats);
    
    %%% 
    for k=1:length(SNR) % for each SNR
        VarNoiseAddQ(k)=VarSigQ/SNR(k);  % variance of added noise to the Q wave
        VarNoiseAddT(k)=VarSigT/SNR(k);  % variance of added noise to the T wave
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
                %%% for ML, fit gaussians on each noisy Q-segs
                MlParamsQ(:,jj,m,k)=GausFit(tt,SegNQ,p0,eval(Qplb),eval(Qpub),optionsMl);
                %%% for Bayesian, fit gaussians on each noisy Q-segs
                BysParamsQ(:,jj,m,k)=GausFit(tt,SegNQ,p0,eval(Qplb),eval(Qpub),PrMeanQ,PrCovQ,VarSigQ/SNR(k),optionsBys);
                
                
                SegNT=signalNT(TSegOn(j):TSegOff(j)); % SOI
                tt=tm(TSegOn(j):TSegOff(j))-Rt(j,2); % SOI time stamp
                p0=[signal(T(j,2)) (Tt(j,3)-Tt(j,1))/5 Tt(j,2)]; % initial points
                %%% for ML, fit gaussians on each noisy T-segs
                MlParamsT(:,jj,m,k)=GausFit(tt,SegNT,p0,eval(Tplb),eval(Tpub),optionsMl);
                %%% for Bayesian, fit gaussians on each noisy T-segs
                BysParamsT(:,jj,m,k)=GausFit(tt,SegNT,p0,eval(Tplb),eval(Tpub),PrMeanT,PrCovT,VarSigT/SNR(k),optionsBys);
            end
            MlErQ(:,:,m,k)=MlParamsQ(:,:,m,k)-PrParamsQ(JJ,:)'; % ML Error from prior
            BysErQ(:,:,m,k)=BysParamsQ(:,:,m,k)-PrParamsQ(JJ,:)'; % Bys Error from prior
            MlErT(:,:,m,k)=MlParamsT(:,:,m,k)-PrParamsT(JJ,:)'; % ML Error from prior
            BysErT(:,:,m,k)=BysParamsT(:,:,m,k)-PrParamsT(JJ,:)'; % Bys Error from prior         
        end
        MlCovQ(:,:,k)=cov(reshape(MlErQ(:,:,:,k),3,[])'); % ML Error covariance matrix
        DetMlCovQ(k)=det(MlCovQ(:,:,k)); % determinant of ML Error covariance matrix
        BysCovQ(:,:,k)=cov(reshape(BysErQ(:,:,:,k),3,[])'); % Bys Error covariance matrix
        DetBysCovQ(k)=det(BysCovQ(:,:,k)); % determinant of Bys Error covariance matrix
        
        MlCovT(:,:,k)=cov(reshape(MlErT(:,:,:,k),3,[])'); % ML Error covariance matrix
        DetMlCovT(k)=det(MlCovT(:,:,k)); % determinant of ML Error covariance matrix
        BysCovT(:,:,k)=cov(reshape(BysErT(:,:,:,k),3,[])'); % Bys Error covariance matrix
        DetBysCovT(k)=det(BysCovT(:,:,k)); % determinant of Bys Error covariance matrix
    end
    close(h)

%% CRLB calculation
%%% pre-allocation
%     MlCrlbApprQ=zeros(3, 3, NumAllBeats, length(SNRdB));
    MlFIApprQ=zeros(3, 3, NumAllBeats, length(SNRdB));
    BysCrlbApprQ=zeros(3, 3, length(SNRdB));
%     MlCrlbApprT=zeros(3, 3, NumAllBeats, length(SNRdB));
    MlFIApprT=zeros(3, 3, NumAllBeats, length(SNRdB));
    BysCrlbApprT=zeros(3, 3, length(SNRdB));
    MlCrlbApprQ_AvrgFI=zeros(3, 3, length(SNRdB));
    MlCrlbApprT_AvrgdFI=zeros(3, 3, length(SNRdB));
    DetMlCrlbApprQ=zeros(length(SNRdB),1);
    DetBysCrlbApprQ=zeros(length(SNRdB),1);
    DetMlCrlbApprT=zeros(length(SNRdB),1);
    DetBysCrlbApprT=zeros(length(SNRdB),1);
    
%     MlCrlbNumQ=zeros(3, 3, NumAllBeats, length(SNRdB));
    MlFINumQ=zeros(3, 3, NumAllBeats, length(SNRdB));
    BysCrlbNumQ=zeros(3, 3, length(SNRdB));
%     MlCrlbNumT=zeros(3, 3, NumAllBeats, length(SNRdB));
    MlFINumT=zeros(3, 3, NumAllBeats, length(SNRdB));
    BysCrlbNumT=zeros(3, 3, length(SNRdB));
    MlCrlbNumQ_AvrgFI=zeros(3, 3, length(SNRdB));
    MlCrlbNumT_AvrgdFI=zeros(3, 3, length(SNRdB));
    DetMlCrlbNumQ=zeros(length(SNRdB),1);
    DetBysCrlbNumQ=zeros(length(SNRdB),1);
    DetMlCrlbNumT=zeros(length(SNRdB),1);
    DetBysCrlbNumT=zeros(length(SNRdB),1);
    
    
    for k=1:length(SNRdB)
        for j=1:NumAllBeats
            [~, MlFINumQ(:,:,j,k)]=GaussCrlbNumeric(tm(QSegOn(j):QSegOff(j))-Rt(j,2), PrParamsQ(j,:), VarNoiseAddQ(k), 1 );
            [~, MlFINumT(:,:,j,k)]=GaussCrlbNumeric(tm(TSegOn(j):TSegOff(j))-Rt(j,2), PrParamsT(j,:), VarNoiseAddT(k), 1 );
            
            [~, MlFIApprQ(:,:,j,k)]=GaussCRLB(PrParamsQ(j,1), PrParamsQ(j,2), Fs, VarNoiseAddQ(k) );
            [~, MlFIApprT(:,:,j,k)]=GaussCRLB(PrParamsT(j,1), PrParamsT(j,2), Fs, VarNoiseAddT(k) );
%             MlCrlbApprT(:,:,j,k)=2*MlCrlbApprT(:,:,j,k);  
            MlFIApprT(:,:,j,k)=.5*MlFIApprT(:,:,j,k); % since we consider half of T wave        
        end
        MlCrlbNumQ_AvrgFI(:,:,k)=inv(mean(MlFINumQ(:,:,:,k),3));
        BysCrlbNumQ(:,:,k)=inv(mean(MlFINumQ(:,:,:,k),3)+inv(PrCovQ));

        MlCrlbNumT_AvrgdFI(:,:,k)=inv(mean(MlFINumT(:,:,:,k),3));
        BysCrlbNumT(:,:,k)=inv(mean(MlFINumT(:,:,:,k),3)+inv(PrCovT));

        DetMlCrlbNumQ(k)=det(MlCrlbNumQ_AvrgFI(:,:,k));
        DetBysCrlbNumQ(k)=det(BysCrlbNumQ(:,:,k));

        DetMlCrlbNumT(k)=det(MlCrlbNumT_AvrgdFI(:,:,k));
        DetBysCrlbNumT(k)=det(BysCrlbNumT(:,:,k));

        MlCrlbApprQ_AvrgFI(:,:,k)=inv(mean(MlFIApprQ(:,:,:,k),3));
        BysCrlbApprQ(:,:,k)=inv(mean(MlFIApprQ(:,:,:,k),3)+inv(PrCovQ));

        MlCrlbApprT_AvrgdFI(:,:,k)=inv(mean(MlFIApprT(:,:,:,k),3));
        BysCrlbApprT(:,:,k)=inv(mean(MlFIApprT(:,:,:,k),3)+inv(PrCovT));

        DetMlCrlbApprQ(k)=det(MlCrlbApprQ_AvrgFI(:,:,k));
        DetBysCrlbApprQ(k)=det(BysCrlbApprQ(:,:,k));

        DetMlCrlbApprT(k)=det(MlCrlbApprT_AvrgdFI(:,:,k));
        DetBysCrlbApprT(k)=det(BysCrlbApprT(:,:,k));
    end

%     save(['Backup and Results\TempQtPrmtrs_' files(i).name(1:end-4)  '.mat'], 'PrParamsQ', 'PrParamsT', 'MlParamsQ', 'BysParamsQ', 'MlParamsT', 'BysParamsT', 'VarNoiseAddQ', 'VarNoiseAddT', 'VarSigQ', 'VarSigT', 'SNR', 'MlCovQ', 'DetMlCovQ', 'BysCovQ', 'DetBysCovQ', 'MlCovT', 'DetMlCovT', 'BysCovT', 'DetBysCovT', 'MlMseQ', 'BysMseQ', 'MlMseT', 'BysMseT', 'MlCrlbQ', 'BysCrlbQ', 'MlCrlbT', 'BysCrlbT', 'MlCrlbQ_AvrgdFI', 'MlCrlbT_AvrgdFI',  'MlFIQ', 'MlFIT', 'DetMlCrlbQ', 'DetBysCrlbQ', 'DetMlCrlbT', 'DetBysCrlbT')    
    save(['Backup and Results\TempQtPrmtrs_' files(i).name(1:end-4)  '.mat']);
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
