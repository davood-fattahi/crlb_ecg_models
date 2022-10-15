clc
clear
close all


%% the constants
N=[]; % the end sample of the signal loaded by rdsamp
N0=1; % the start sample of the signal loaded by rdsamp
Chnl = 1; % ECG channel 
HrTr= 70/60; % Hz, threshold for heart rate, used in r peak detector
SNRdB=(-20:5:40); % vector of SNRs in decibel
SNR=10.^(SNRdB/10); % SNRs
NumGaus=1; % number of Gaussians utilized in modeling ecg segments
w1=.25; % window length of the first median filter used in base line removing
w2=.30; % window length of the second median filter used in base line removing
NumBeats=[]; % number of the beats envolved in noisy signals' parameter estimation stage
NumRuns=3; % number of repeats in noisy signals' parameter estimation stage
AllBeats=[]; % number of the beats envolved in prior parameter estimation stage (clean signals)

%%% option structure for non linear least squre optimization, respectively
%%% in ML estimation and Bayesian estimation:
options = struct('SpecifyObjectiveGradient',true);
% options = struct('SpecifyObjectiveGradient',true, 'FunctionTolerance', 1e-8, 'Display','final-detailed' ,'MaxFunctionEvaluations',10000, 'MaxIterations', 10000, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12);
% options=struct('SpecifyObjectiveGradient',true);


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

AllPrParamsQ=[];
AllPrParamsT=[];
ii=[];
for i= 1:NumCases
    %% load the record
    oldFolder=cd('..\mcode'); % go to the wfdb toolbox root
    AdrsNm=[DA  '\' files(i).name(1:end-4)]; % address of the dataset
    [signal,Fs,tm]=rdsamp(AdrsNm,[],N,N0); % loading the data
    
    
    % Import the annotation
    try
    [annatr,anntypeatr,subtypeatr,chanatr,numatr,commentsatr]=rdann(AdrsNm, 'atr', [], N, N0);
    [annpu0,anntypepu0,subtypepu0,chanpu0,numpu0,commentspu0]=rdann(AdrsNm, 'pu0', [], N, N0);
    catch ME
        warning(['Case No. ' num2str(i) ', ' files(i).name(1:end-4) ' is skiped. See the error bellow:' newline ME.message]);
        cd(oldFolder)
        continue
    end
    cd(oldFolder) % go back to the main folder
    
    
    %% finding the beats with a normal QRS and a normal T 
    TNormalT=annpu0(anntypepu0=='t'&numpu0==0);
    RNormalQRS=annatr(anntypeatr=='N');
    A=(repmat(TNormalT,1,length(annatr))-annatr')<0;
    [~,c]=find(([A(:,2:end) true(size(A,1),1)]-A)==1);
    RNormalT=annatr(c);
    RNormalQRST=RNormalT(ismember(RNormalT,RNormalQRS));
    TNormalQRST=TNormalT(ismember(RNormalT,RNormalQRS));
    
    if isempty(RNormalQRST)
        continue
    end
    %% baseline wandering removal,
    signal=signal(:,Chnl);    
    signal=signal-(BaseLine1(BaseLine1(signal', round(w1*Fs), 'md'), round(w2*Fs), 'mn'))';
  
    %% adjusting the AllBeats by the detected R peaks:
    Rp=RNormalQRST;
    if isempty(AllBeats) || (AllBeats >length(Rp))  %#ok<BDSCI>
        NumAllBeats=length(Rp);
    else
        NumAllBeats=AllBeats;
    end
    
    AllBeat(
    
    %% extract approx. SOI
    [Q, T] =QTfidapprox(signal, Rp, Fs);
    Rpt=tm(Rp); Qt=tm(Q)-Rpt; Tt=tm(T)-Rpt; % getting time stamps of the fid points
    
    %%% adjusting the start and the end of SOI
    QSegOn=Q(:,Qsegon(1))+floor(Qsegon(2)*Fs/1000); 
    QSegOff=Q(:,Qsegoff(1))+floor(Qsegoff(2)*Fs/1000); 
    TSegOn=T(:,Tsegon(1))+floor(Tsegon(2)*Fs/1000);
    TSegOff=T(:,Tsegoff(1))+floor(Tsegoff(2)*Fs/1000);
    
%     % plot the fid points on the signal
%     figure
%     plot(tm,signal); hold on
%     plot(tm(Q(~isnan(Q))),signal(Q(~isnan(Q))),'*')
%     plot(tm(Rp(~isnan(Rp))),signal(Rp(~isnan(Rp))),'*')
%     plot(tm(T(~isnan(T))),signal(T(~isnan(T))),'*')    
% %     saveas(gcf,'ApprFP.fig')


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
        tt=tm(QSegOn(j):QSegOff(j))-Rpt(j); % SOI time stamp
        
        % fit gssns on each clean Q-wave
        p0=[signal(Q(j,2)) (Qt(j,3)-Qt(j,1))/5 Qt(j,2)]; % initial condition
        PrParamsQ(j,:)=GausFit(tt,SegQ,p0,eval(Qplb),eval(Qpub),options);
        
%         % polt the evaluated gaussians on the signal
%         SegEvQ=GausVal(tt,PrParamsQ(j,:));
%         hold on;
%         plot(tt+Rpt(j),SegEvQ,'r-')
%        
        %% T wave model fitting      
        SegT=signal(TSegOn(j):TSegOff(j)); % Segment Of Intrest (SOI)
        SegsT(LT+1:LT+length(SegT))=SegT(:); % All the SOIs
        LT=LT+length(SegT); % All the SOIs length
        tt=tm(TSegOn(j):TSegOff(j))-Rpt(j); % SOI time stamp
        
        % fit gssns on each clean T-halfwave
        p0=[signal(T(j,2)) (Tt(j,3)-Tt(j,1))/5 Tt(j,2)]; % initial condition
        PrParamsT(j,:)=GausFit(tt,SegT,p0,eval(Tplb),eval(Tpub),options);
        
%         % polt the evaluated gaussians on the signal
%         tt=tm(T(j,1):T(j,3))-Rpt(j);
%         SegEvT=GausVal(tt(:),PrParamsT(j,:));
%         hold on
%         plot(tt(:)+Rpt(j),SegEvT,'r-')

    end
    close(h)
    VarSigQ=var(SegsQ(:)); % variance of the Q signal
    VarSigT=var(SegsT(:)); % variance of the T signal
    AllPrParamsQ(end+1:end+NumAllBeats,:)=PrParamsQ;
    AllPrParamsT(end+1:end+NumAllBeats,:)=PrParamsT;
    save(['Backup and Results\PriorQtPrmtrs_' files(i).name(1:end-4)  '.mat']);
end

    %% estimate the prior for the params
    PrMeanQ=mean(AllPrParamsQ,1);
    PrCovQ=cov(AllPrParamsQ);
    PrMeanT=mean(AllPrParamsT,1);
    PrCovT=cov(AllPrParamsT);
    
    %% Stage 2- estimating the parameters from noisy signals
    h = waitbar(0,['Case No. ' num2str(i) ', Estimating parameters for noisy signals, please wait ...']);
for i=ii
    load(['Backup and Results\PriorQtPrmtrs_' files(i).name(1:end-4)  '.mat']);

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
                tt=tm(QSegOn(j):QSegOff(j))-Rpt(j); % SOI time stamp
                p0=[signal(Q(j,2)) (Qt(j,3)-Qt(j,1))/5 Qt(j,2)]; % initial points
                %%% for ML, fit gaussians on each noisy Q-segs
                MlParamsQ(:,jj,m,k)=GausFit(tt,SegNQ,p0,eval(Qplb),eval(Qpub),options);
                %%% for Bayesian, fit gaussians on each noisy Q-segs
                BysParamsQ(:,jj,m,k)=GausFit(tt,SegNQ,p0,eval(Qplb),eval(Qpub),PrMeanQ,PrCovQ,VarSigQ/SNR(k),options);
                
                
                SegNT=signalNT(TSegOn(j):TSegOff(j)); % SOI
                tt=tm(TSegOn(j):TSegOff(j))-Rpt(j); % SOI time stamp
                p0=[signal(T(j,2)) (Tt(j,3)-Tt(j,1))/5 Tt(j,2)]; % initial points
                %%% for ML, fit gaussians on each noisy T-segs
                MlParamsT(:,jj,m,k)=GausFit(tt,SegNT,p0,eval(Tplb),eval(Tpub),options);
                %%% for Bayesian, fit gaussians on each noisy T-segs
                BysParamsT(:,jj,m,k)=GausFit(tt,SegNT,p0,eval(Tplb),eval(Tpub),PrMeanT,PrCovT,VarSigT/SNR(k),options);
            end
            MlErQ(:,:,m,k)=MlParamsQ(:,:,m,k)-PrParamsQ(JJ,:)'; % ML Error from prior
            BysErQ(:,:,m,k)=BysParamsQ(:,:,m,k)-PrParamsQ(JJ,:)'; % Bys Error from prior
            MlErT(:,:,m,k)=MlParamsT(:,:,m,k)-PrParamsT(JJ,:)'; % ML Error from prior
            BysErT(:,:,m,k)=BysParamsT(:,:,m,k)-PrParamsT(JJ,:)'; % Bys Error from prior         
        end
    end
    close(h)

%% CRLB calculation
%%% pre-allocation
    MlFIApprQ=zeros(3, 3, NumAllBeats, length(SNRdB));
    MlFIApprT=zeros(3, 3, NumAllBeats, length(SNRdB));
   
    MlFINumQ=zeros(3, 3, NumAllBeats, length(SNRdB));
    MlFINumT=zeros(3, 3, NumAllBeats, length(SNRdB));
  
    
    for k=1:length(SNRdB)
        for j=1:NumAllBeats
            [~, MlFINumQ(:,:,j,i,k)]=GaussCrlbNumeric(tm(QSegOn(j):QSegOff(j))-Rpt(j), PrParamsQ(j,:), VarNoiseAddQ(k), 1 );
            [~, MlFINumT(:,:,j,i,k)]=GaussCrlbNumeric(tm(TSegOn(j):TSegOff(j))-Rpt(j), PrParamsT(j,:), VarNoiseAddT(k), 1 );
            
            [~, MlFIApprQ(:,:,j,i,k)]=GaussCRLB(PrParamsQ(j,1), PrParamsQ(j,2), Fs, VarNoiseAddQ(k) );
            [~, MlFIApprT(:,:,j,i,k)]=GaussCRLB(PrParamsT(j,1), PrParamsT(j,2), Fs, VarNoiseAddT(k) );
%             MlCrlbApprT(:,:,j,k)=2*MlCrlbApprT(:,:,j,i,k);  
            MlFIApprT(:,:,j,i,k)=.5*MlFIApprT(:,:,j,i,k); % since we consider half of T wave        
        end
    end

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
