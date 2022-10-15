clc
clear
% close all


%% the constants
Chnl = [5 2 5 11 4  4  6 2 1 2  4  4 1 7  3 1 12 1 4 1 1 1  4 6;
        1 1 1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1 1 1 1 -1 1 1 1 -1 1]; % ECG channels and their polarity.
% phase0 = 0;   % desired phase shift
% bins = 250;   % number of phase bins
SynthEcgFs=1000;
pstd=.03; % parameters relative standard deviation around their mean.
rrstd=.04; % rr interval relative standard deviation around its mean.
SynthSigLength=1*.5*60; % length of synthetic ecg (sec)
HrTr= 70/60; % Hz, threshold for heart rate, used in r peak detector
SNRdB=(-20:20:40); % vector of SNRs in decibel
SNR=10.^(SNRdB/10); % SNRs
NumGaus=1; % number of Gaussians utilized in modeling ecg segments
% w1=.25; % window length of the first median filter used in base line removing
% w2=.30; % window length of the second median filter used in base line removing
PostNumBeats=[]; % max number of the beats envolved in noisy signals' parameter estimation stage
NumRuns=1; % number of repeats in noisy signals' parameter estimation stage
PriNumBeats=[]; % max number of the beats envolved in prior parameter estimation stage (clean signals)

%%% option structure for non linear least squre optimization, respectively
%%% in ML estimation and Bayesian estimation:
options = struct('SpecifyObjectiveGradient',true);
% options = struct('SpecifyObjectiveGradient',true, 'FunctionTolerance', 1e-10, 'Display','final-detailed' ,'MaxFunctionEvaluations',10000, 'MaxIterations', 10000, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12);
% options=struct('SpecifyObjectiveGradient',true);

Qplb='[-5 0 tm(QSegOn(j))  -10 0 tm(QSegOff(j))]'; % optimization conditions: lower bound for the Q wave parameters
Qpub='[5  3*(tm(QSegOff(j))-tm(QSegOn(j))) tm(QSegOff(j)) 10 3*(tm(RSegOff(j))-tm(QSegOff(j))) tm(RSegOff(j))]'; % optimization conditions: upper bound for the Q wave parameters

Tplb='[-5 0 tm(TSegOn(j))]';% optimization conditions: lower bound for the T wave parameters
Tpub='[5 3*(tm(TSegOff(j))-tm(TSegOn(j))) tm(TSegOff(j))]';% optimization conditions: upper bound for the T wave parameters


%% load the case names in the directory
DA='..\mcode\database\ptb-gaussian-parameters';
files=dir([DA '\*.mat']); 
NumCases=length(files);

ii=1;
AllBeats=[]; Rp=[]; Qp=[]; Tp=[]; CaseNum=[];
QSegOn=[]; QSegOff=[]; RSegOff=[]; TSegOn=[]; TSegOff=[]; 
%% stage 1- gathering the approporiat beats
h=waitbar(0,['gathering the appropriate beats, please wait ...']);
for i=1:NumCases %[1, 3, 4, 5, 6, 13, 17, 19, 20]
    waitbar(i/NumCases);

    %% load the parameters
    GausParam=load([DA '\' files(i).name]);
    pmean=[1e-3*Chnl(2,i).*GausParam.params{Chnl(1,i)}.a; GausParam.params{Chnl(1,i)}.b; GausParam.params{Chnl(1,i)}.theta];
    
    %% resemble the signal
    tm=1/SynthEcgFs:1/SynthEcgFs:SynthSigLength;
    rng(i) % fixing the seed for random variation in beats
    [signal, ~]=ECGResembler(tm, pmean, pstd, 1/HrTr, rrstd);
    signal=signal(:);
    
    %% baseline wandering removal,
%     signal=signal-(BaseLine1(BaseLine1(signal', round(w1*SynthEcgFs), 'md'), round(w2*SynthEcgFs), 'mn'))';
    
    %% extract the R peaks using peak detector
    Rpeaks = PeakDetection(signal(:),HrTr/SynthEcgFs); 
    Rpeaks = find(Rpeaks); NumBeats=length(Rpeaks);
    [~, QQ, RR, ~, TT] = ecgWavesSoI_nrmlBeat(signal, Rpeaks, 1);
    

    Rp(end+1:end+NumBeats-2,:) = RR(2:end-1,2)+floor(.4*SynthEcgFs)+1;
    RSegOff(end+1:end+NumBeats-2,:) = RR(2:end-1,3)+floor(.4*SynthEcgFs)+1;
    
    QSegOn(end+1:end+NumBeats-2,:) = QQ(2:end-1,1)+floor(.4*SynthEcgFs)+1;
    Qp(end+1:end+NumBeats-2,:) = QQ(2:end-1,2)+floor(.4*SynthEcgFs)+1;
    QSegOff(end+1:end+NumBeats-2,:) = QQ(2:end-1,3)+floor(.4*SynthEcgFs)+1;
    
    TSegOn(end+1:end+NumBeats-2,:) = TT(2:end-1,2)- fix(0.02*SynthEcgFs) + floor(.4*SynthEcgFs)+1;
    Tp(end+1:end+NumBeats-2,:) = TT(2:end-1,2)+floor(.4*SynthEcgFs)+1;
    TSegOff(end+1:end+NumBeats-2,:) = TT(2:end-1,3) + floor(.4*SynthEcgFs)+1;
    
    CaseNum(end+1:end+NumBeats-2,:)=i;
    
    b=-floor(.4*SynthEcgFs):floor(.6*SynthEcgFs);
    bb = Rpeaks' + b(ones(size(Rpeaks))',:);
    AllBeats(end+1:end+NumBeats-2,:)=reshape(signal(bb(2:end-1,:)),[],length(b));
    
    % Note: Q, Rp and T may be the same for all the beats, but we allocate
    % them separately for each beat for possible changes in the future. 

end
close(h);
NumAllBeats=length(Rp);

%% adjusting the PriNumBeats by NumAllBeats:
if isempty(PriNumBeats) || (PriNumBeats > NumAllBeats)  %#ok<BDSCI>
    PriNumBeats=NumAllBeats;
end

%%% selecting beats randomly
rng(2) % fixing the seed for random permuting
JJ=randperm(NumAllBeats,PriNumBeats);
AllBeats=AllBeats(JJ,:); 

Rp=Rp(JJ);
RSegOff = RSegOff(JJ);

QSegOn = QSegOn(JJ);
Qp = Qp(JJ);
QSegOff = QSegOff(JJ);

TSegOn = TSegOn(JJ);
TSegOff = TSegOff(JJ);
Tp = Tp(JJ);
CaseNum=CaseNum(JJ);

save('Backup and Results\SynthDataQRT_InitialPool.mat');

%% stage 2- computing the prior

% Pre-allocating ...
PrParamsQ=zeros(PriNumBeats,6);
PrParamsT=zeros(PriNumBeats,3);
SegsQ=nan(size(AllBeats));
SegsT=nan(size(AllBeats));

% wait bar
h = waitbar(0,'Estimating parameters for clean signals, please wait ...');
for j=1:PriNumBeats  % for each beat ...
    waitbar(j/PriNumBeats)
    % extract approx. SOI
    tm=((1:size(AllBeats,2))-Rp(j))/SynthEcgFs; % time stamp; the R peak is the reference time (zero).        

    % fit gssns on each clean Q-wave
    SegQR=AllBeats(j,QSegOn(j):RSegOff(j)); % Segment Of Intrest (SOI)
    SegsQ(j,QSegOn(j):QSegOff(j))=AllBeats(j,QSegOn(j):QSegOff(j)); % All the SOIs
    p0=[AllBeats(j,Qp(j)) (tm(RSegOff(j))-tm(QSegOn(j)))/5 tm(Qp(j)) ...
        AllBeats(j,Rp(j)) (tm(RSegOff(j))-tm(QSegOn(j)))/5 tm(Rp(j)) ]; % initial condition
    p=GausFit(tm(QSegOn(j):RSegOff(j)),SegQR,p0,eval(Qplb),eval(Qpub),options);
    [~,I]=sort(p([3 6]));
    PrParamsQ(j,:) = [p((I(1)-1)*3+1:(I(1)-1)*3+3); p((I(2)-1)*3+1:(I(2)-1)*3+3)] ;
    
    
    % fit gssns on each clean T-halfwave
    SegT=AllBeats(j,TSegOn(j):TSegOff(j)); % Segment Of Intrest (SOI)
    SegsT(j,TSegOn(j):TSegOff(j))=SegT(:); % All the SOIs
    p0=[AllBeats(j,Tp(j)) (tm(TSegOff(j))-tm(TSegOn(j)))/5 tm(Tp(j))]; % initial condition
    PrParamsT(j,:)=GausFit(tm(TSegOn(j):TSegOff(j)),SegT,p0,eval(Tplb),eval(Tpub),options);

%     % polt the evaluated gaussians on the signal
%     hold off; plot(tm,AllBeats(j,:)); hold on
%     plot(tm(QSegOn(j):RSegOff(j)),GausVal(tm(QSegOn(j):RSegOff(j)),PrParamsQ(j,:)),'r-')    
%     plot(tm(TSegOn(j):TSegOff(j)),GausVal(tm(TSegOn(j):TSegOff(j)),PrParamsT(j,:)),'r-')
end
close(h)

% estimate the prior for the params
PrMeanQ=mean(PrParamsQ,1);
PrCovQ=cov(PrParamsQ);
PrMeanT=mean(PrParamsT,1);
PrCovT=cov(PrParamsT);
VarSigQ=nan(PriNumBeats,1);
VarSigT=nan(PriNumBeats,1);

for i = unique(CaseNum)'
    VarSigQ(CaseNum==i)=var(SegsQ(CaseNum==i,:),[],'all','omitnan');
    VarSigT(CaseNum==i)=var(SegsT(CaseNum==i,:),[],'all','omitnan');
end
% VarSigQ(:)=var(SegsQ,[],'all','omitnan');
% VarSigT(:)=var(SegsT,[],'all','omitnan');

save('Backup and Results\SynthDataQRT_PriorQtPrmtrs.mat');

%% Stage 3- estimating the parameters from noisy signals
h = waitbar(0,'Estimating parameters for noisy signals, please wait ...');
%%% adjusting the number of noisy beats
if isempty(PostNumBeats)
    PostNumBeats=PriNumBeats;
end
if PostNumBeats > PriNumBeats
    PostNumBeats = PriNumBeats;
end

%%% pre allocating ...
MlParamsQ=zeros(2*3, PostNumBeats, NumRuns, length(SNR));
MlErQ=zeros(3, PostNumBeats, NumRuns, length(SNR));   

BysParamsQ=zeros(2*3, PostNumBeats, NumRuns, length(SNR));
BysErQ=zeros(3, PostNumBeats, NumRuns, length(SNR));

MlParamsT=zeros(1*3, PostNumBeats, NumRuns, length(SNR)); 
MlErT=zeros(size(MlParamsT)); 

BysParamsT=zeros(1*3, PostNumBeats, NumRuns, length(SNR));
BysErT=zeros(size(BysParamsT));     

MlCovQ=zeros(3,3,length(SNR));
BysCovQ=zeros(3,3,length(SNR));

DetMlCovQ=zeros(size(SNR)); 
DetBysCovQ=zeros(size(SNR)); 

MlCovT=zeros(3,3,length(SNR));
BysCovT=zeros(3,3,length(SNR));

DetMlCovT=zeros(size(SNR)); 
DetBysCovT=zeros(size(SNR)); 
 
%%% select the noisy beats randomly
rng(3) % fixing the seed for random permuting
JJ=randperm(PriNumBeats,PostNumBeats);

%%% 
for k=1:length(SNR) % for each SNR
    for m=1:NumRuns % for each run
        jj=0; % initializing
        for j=JJ
            waitbar(((k-1)*(NumRuns)*length(JJ)+(m-1)*length(JJ)+jj)/(length(SNR)*NumRuns*length(JJ))) % wait bar ...
            jj=jj+1; % number of beats counter
    
            tm=((1:size(AllBeats,2))-Rp(j))/SynthEcgFs; % time stamp; the R peak is the reference time (zero).        
            rng(k+m+j); % fix the seed for noise generating
            SegNQR=AllBeats(j,QSegOn(j):RSegOff(j))+randn(1,RSegOff(j)-QSegOn(j)+1).*sqrt(VarSigQ(j)/SNR(k)); % adding noise to the signal;
            tt=tm(QSegOn(j):RSegOff(j)); % SOI time stamp
            p0=[AllBeats(j,Qp(j)) (tm(RSegOff(j))-tm(QSegOn(j)))/5 tm(Qp(j)) ...
                AllBeats(j,Rp(j)) (tm(RSegOff(j))-tm(QSegOn(j)))/5 tm(Rp(j)) ]; % initial condition
            %%% for ML, fit gaussians on each noisy Q-segs
            p=GausFit(tt,SegNQR,p0,eval(Qplb),eval(Qpub),options);
            [~,I]=sort(p([3 6])); MlParamsQ(:,jj,m,k)=[p((I(1)-1)*3+1:(I(1)-1)*3+3); p((I(2)-1)*3+1:(I(2)-1)*3+3)] ;
            %%% for Bayesian, fit gaussians on each noisy Q-segs
            p=GausFit(tt,SegNQR,p0,eval(Qplb),eval(Qpub),PrMeanQ,PrCovQ,VarSigQ(j)/SNR(k),options);
            [~,I]=sort(p([3 6])); BysParamsQ(:,jj,m,k)=[p((I(1)-1)*3+1:(I(1)-1)*3+3); p((I(2)-1)*3+1:(I(2)-1)*3+3)] ;

            rng(k+m+j+1); % fix the seed for noise generating
            SegNT=AllBeats(j,TSegOn(j):TSegOff(j)) +randn(1,TSegOff(j)-TSegOn(j)+1).*sqrt(VarSigT(j)/SNR(k)); % adding noise to the signal
            tt=tm(TSegOn(j):TSegOff(j)); % SOI time stamp
            p0=[max(SegNT) length(SegNT)/3/SynthEcgFs tm(Tp(j))]; % initial points
            %%% for ML, fit gaussians on each noisy T-segs
            MlParamsT(:,jj,m,k)=GausFit(tt,SegNT,p0,eval(Tplb),eval(Tpub),options);
            %%% for Bayesian, fit gaussians on each noisy T-segs
            BysParamsT(:,jj,m,k)=GausFit(tt,SegNT,p0,eval(Tplb),eval(Tpub),PrMeanT,PrCovT,VarSigT(j)/SNR(k),options);

            
%             % polt the evaluated gaussians on the signal
%             hold off; plot(tm,AllBeats(j,:)); hold on
%             plot(tm(QSegOn(j):RSegOff(j)),GausVal(tm(QSegOn(j):RSegOff(j)),MlParamsQ(:,jj,m,k)),'r-')    
%             plot(tm(TSegOn(j):TSegOff(j)),GausVal(tm(TSegOn(j):TSegOff(j)),MlParamsT(:,jj,m,k)),'r-') 
        end
        MlErQ(:,:,m,k)=MlParamsQ(1:3,:,m,k)-PrParamsQ(JJ,1:3)'; % ML Error from prior
        BysErQ(:,:,m,k)=BysParamsQ(1:3,:,m,k)-PrParamsQ(JJ,1:3)'; % Bys Error from prior
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
%     MlCrlbApprQ=zeros(3, 3, PriNumBeats, length(SNRdB));
% MlFIApprQ=zeros(6, 6, PriNumBeats, length(SNRdB));
% BysCrlbApprQ=zeros(6, 6, length(SNRdB));
% %     MlCrlbApprT=zeros(3, 3, PriNumBeats, length(SNRdB));
% MlFIApprT=zeros(3, 3, PriNumBeats, length(SNRdB));
% BysCrlbApprT=zeros(3, 3, length(SNRdB));
% MlCrlbApprQ_AvrgdFI=zeros(6, 6, length(SNRdB));
% MlCrlbApprT_AvrgdFI=zeros(3, 3, length(SNRdB));
% DetMlCrlbApprQ=zeros(length(SNRdB),1);
% DetBysCrlbApprQ=zeros(length(SNRdB),1);
% DetMlCrlbApprT=zeros(length(SNRdB),1);
% DetBysCrlbApprT=zeros(length(SNRdB),1);

%     MlCrlbNumQ=zeros(3, 3, PriNumBeats, length(SNRdB));
MlFINumQ=zeros(6, 6, PriNumBeats, length(SNRdB));
BysCrlbNumQ=zeros(3, 3, length(SNRdB));
%     MlCrlbNumT=zeros(3, 3, PriNumBeats, length(SNRdB));
MlFINumT=zeros(3, 3, PriNumBeats, length(SNRdB));
BysCrlbNumT=zeros(3, 3, length(SNRdB));
MlCrlbNumQ_AvrgdFI=zeros(3, 3, length(SNRdB));
MlCrlbNumT_AvrgdFI=zeros(3, 3, length(SNRdB));
DetMlCrlbNumQ=zeros(length(SNRdB),1);
DetBysCrlbNumQ=zeros(length(SNRdB),1);
DetMlCrlbNumT=zeros(length(SNRdB),1);
DetBysCrlbNumT=zeros(length(SNRdB),1);


for k=1:length(SNRdB)
    for j=1:PriNumBeats
        [~, MlFINumQ(:,:,j,k)]=GaussCrlbNumeric(tm(QSegOn(j):QSegOff(j)), PrParamsQ(j,:), VarSigQ(j)/SNR(k), 1 );
        [~, MlFINumT(:,:,j,k)]=GaussCrlbNumeric(tm(TSegOn(j):TSegOff(j)), PrParamsT(j,:), VarSigT(j)/SNR(k), 1 );

%         [~, MlFIApprQ(:,:,j,k)]=GaussCRLB(PrParamsQ(j,1), PrParamsQ(j,2), SynthEcgFs, VarSigQ(j)/SNR(k) );
%         [~, MlFIApprT(:,:,j,k)]=GaussCRLB(PrParamsT(j,1), PrParamsT(j,2), SynthEcgFs, VarSigT(j)/SNR(k) );
%             MlCrlbApprT(:,:,j,k)=2*MlCrlbApprT(:,:,j,k);  
%         MlFIApprT(:,:,j,k)=.5*MlFIApprT(:,:,j,k); % since we consider half of T wave        
    end
    MlCrlbNumQ_AvrgdFI(:,:,k)=inv(mean(MlFINumQ(1:3,1:3,:,k),3));
    BysCrlbNumQ(:,:,k)=inv(mean(MlFINumQ(1:3,1:3,:,k),3)+inv(PrCovQ(1:3,1:3)));

    DetMlCrlbNumQ(k)=det(MlCrlbNumQ_AvrgdFI(1:3,1:3,k));
    DetBysCrlbNumQ(k)=det(BysCrlbNumQ(1:3,1:3,k));

%     MlCrlbApprQ_AvrgdFI(:,:,k)=inv(mean(MlFIApprQ(:,:,:,k),3));
%     BysCrlbApprQ(:,:,k)=inv(mean(MlFIApprQ(:,:,:,k),3)+inv(PrCovQ));

%     DetMlCrlbApprQ(k)=det(MlCrlbApprQ_AvrgdFI(:,:,k));
%     DetBysCrlbApprQ(k)=det(BysCrlbApprQ(:,:,k));
    
    MlCrlbNumT_AvrgdFI(:,:,k)=inv(mean(MlFINumT(:,:,:,k),3));
    BysCrlbNumT(:,:,k)=inv(mean(MlFINumT(:,:,:,k),3)+inv(PrCovT));

    DetMlCrlbNumT(k)=det(MlCrlbNumT_AvrgdFI(:,:,k));
    DetBysCrlbNumT(k)=det(BysCrlbNumT(:,:,k));


%     MlCrlbApprT_AvrgdFI(:,:,k)=inv(mean(MlFIApprT(:,:,:,k),3));
%     BysCrlbApprT(:,:,k)=inv(mean(MlFIApprT(:,:,:,k),3)+inv(PrCovT));

%     DetMlCrlbApprT(k)=det(MlCrlbApprT_AvrgdFI(:,:,k));
%     DetBysCrlbApprT(k)=det(BysCrlbApprT(:,:,k));
end
save('Backup and Results\SynthDataQrT_Results.mat');

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
% semilogy(SNRdB,DetMlCovQ,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
% semilogy(SNRdB,DetMlCrlbNumQ,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
% semilogy(SNRdB,DetBysCovQ,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
% semilogy(SNRdB,DetBysCrlbNumQ,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
% xlabel 'SNR (dB)'; ylabel 'Determinant (mV^2.Sec^4)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% % title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
% saveas(gcf,'Backup and Results\DetQ_SynthData.fig')
% saveas(gcf,'Backup and Results\DetQ_SynthData.eps','epsc')
% 
% 
% figure
% semilogy(SNRdB,DetMlCovT,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
% semilogy(SNRdB,DetMlCrlbNumT,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
% semilogy(SNRdB,DetBysCovT,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
% semilogy(SNRdB,DetBysCrlbNumT,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
% xlabel 'SNR (dB)'; ylabel 'Determinant (mV^2.Sec^4)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% % title 'Determinants of Error Cov. and CRLB matrices for T wave Gaus. parameters'
% saveas(gcf,'Backup and Results\DetT_SynthData.fig')
% saveas(gcf,'Backup and Results\DetT_SynthData.eps','epsc')


%% QT interval error and CRLB
% clear
% load('Backup and Results\results heavy run\SynthDataQT_Results.mat');
% close all
beta=3;
PrToff=PrParamsT(:,3) + beta*PrParamsT(:,2);
PrQon=PrParamsQ(:,3) - beta*PrParamsQ(:,2);
PrQtInt=PrToff - PrQon;


MlParamsToff=reshape(MlParamsT(3,:,:,:) + beta*MlParamsT(2,:,:,:),length(JJ),NumRuns,length(SNR));
MlParamsQon=reshape(MlParamsQ(3,:,:,:) - beta*MlParamsQ(2,:,:,:),length(JJ),NumRuns,length(SNR));
MlParamsQtInt=MlParamsToff - MlParamsQon;

BysParamsToff=reshape(BysParamsT(3,:,:,:) + beta*BysParamsT(2,:,:,:),length(JJ),NumRuns,length(SNR));
BysParamsQon=reshape(BysParamsQ(3,:,:,:) - beta*BysParamsQ(2,:,:,:),length(JJ),NumRuns,length(SNR));
BysParamsQtInt=BysParamsToff-BysParamsQon;


MlErQon=MlParamsQon-PrQon(JJ); % ML Error 
BysErQon=BysParamsQon-PrQon(JJ); % Bys Error 

MlErToff=MlParamsToff-PrToff(JJ); % ML Error 
BysErToff=BysParamsToff-PrToff(JJ); % Bys Error 

MlErQtInt=MlParamsQtInt-PrQtInt(JJ); % ML Error 
BysErQtInt=BysParamsQtInt-PrQtInt(JJ); % Bys Error 


MlCovQon=zeros(size(SNR)); BysCovQon=zeros(size(SNR)); MlCrlbQon=zeros(size(SNR)); BysCrlbQon=zeros(size(SNR));
for k=1:length(SNR)
    MlCovQon(k)=var(MlErQon(:,:,k),[],'all'); % ML Error variance
    BysCovQon(k)=var(BysErQon(:,:,k),[],'all'); % Bys Error varinace
    MlCrlbQon(k)=[0 -beta 1]*MlCrlbNumQ_AvrgdFI(1:3,1:3,k)*[0 -beta 1]';
    BysCrlbQon(k)=[0 -beta 1]*BysCrlbNumQ(1:3,1:3,k)*[0 -beta 1]';
end


MlCovToff=zeros(size(SNR)); BysCovToff=zeros(size(SNR)); MlCrlbToff=zeros(size(SNR)); BysCrlbToff=zeros(size(SNR));
for k=1:length(SNR)
    MlCovToff(k)=var(MlErToff(:,:,k),[],'all'); % ML Error variance
    BysCovToff(k)=var(BysErToff(:,:,k),[],'all'); % Bys Error varinace
    MlCrlbToff(k)= [0 beta 1]*MlCrlbNumT_AvrgdFI(:,:,k)*[0 beta 1]';
    BysCrlbToff(k)=  [0 beta 1]*BysCrlbNumT(:,:,k)*[0 beta 1]';
end


MlCovQtInt=zeros(size(SNR)); BysCovQtInt=zeros(size(SNR)); MlCrlbQtInt=zeros(size(SNR)); BysCrlbQtInt=zeros(size(SNR));
for k=1:length(SNR)
    MlCovQtInt(k)=var(MlErQtInt(:,:,k),[],'all'); % ML Error variance
    BysCovQtInt(k)=var(BysErQtInt(:,:,k),[],'all'); % Bys Error varinace
    MlCrlbQtInt(k)=[0 -beta 1]*MlCrlbNumQ_AvrgdFI(1:3,1:3,k)*[0 -beta 1]' + [0 beta 1]*MlCrlbNumT_AvrgdFI(1:3,1:3,k)*[0 beta 1]';
    BysCrlbQtInt(k)=[0 -beta 1]*BysCrlbNumQ(1:3,1:3,k)*[0 -beta 1]' + [0 beta 1]*BysCrlbNumT(1:3,1:3,k)*[0 beta 1]';
end
% 
% 
% figure
% semilogy(SNRdB,MlCovQtInt,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
% semilogy(SNRdB,MlCrlbQtInt,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
% semilogy(SNRdB,BysCovQtInt,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
% semilogy(SNRdB,BysCrlbQtInt,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
% xlabel 'SNR (dB)'; ylabel 'Error varinace (Sec^2)'; grid on
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% % title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
% saveas(gcf,'Backup and Results\QTintVarEr_SynthData.fig')
% saveas(gcf,'Backup and Results\QTintVarEr_SynthData.eps','epsc')
% 
% 
% 
% figure
% semilogy(SNRdB,MlCovQon,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
% semilogy(SNRdB,MlCrlbQon,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
% semilogy(SNRdB,BysCovQon,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
% semilogy(SNRdB,BysCrlbQon,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
% xlabel 'SNR (dB)'; ylabel 'Error varinace (Sec^2)'; grid on
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% % title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
% saveas(gcf,'Backup and Results\QonVarEr_SynthData.fig')
% saveas(gcf,'Backup and Results\QonVarEr_SynthData.eps','epsc')
% 
% 
% 
% figure
% semilogy(SNRdB,MlCovToff,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
% semilogy(SNRdB,MlCrlbToff,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
% semilogy(SNRdB,BysCovToff,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
% semilogy(SNRdB,BysCrlbToff,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
% xlabel 'SNR (dB)'; ylabel 'Error varinace (Sec^2)'; grid on
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% % title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
% saveas(gcf,'Backup and Results\ToffVarEr_SynthData.fig')
% saveas(gcf,'Backup and Results\ToffVarEr_SynthData.eps','epsc')
% 


figure
semilogy(SNRdB,MlCovQtInt.^.5,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,MlCrlbQtInt.^.5,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,BysCovQtInt.^.5,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,BysCrlbQtInt.^.5,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Error STD (Sec)'; grid on
legend('ML RMSE','ML CRLB^{1/2}',' BYS RMSE','BYS CRLB^{1/2}')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'Backup and Results\QTintStdEr_SynthData.fig')
saveas(gcf,'Backup and Results\QTintStdEr_SynthData.eps','epsc')



figure
semilogy(SNRdB,MlCovQon.^.5,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,MlCrlbQon.^.5,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,BysCovQon.^.5,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,BysCrlbQon.^.5,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Error STD (Sec)'; grid on
legend('ML RMSE','ML CRLB^{1/2}',' BYS RMSE','BYS CRLB^{1/2}')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'Backup and Results\QonStdEr_SynthData.fig')
saveas(gcf,'Backup and Results\QonStdEr_SynthData.eps','epsc')



figure
semilogy(SNRdB,MlCovToff.^.5,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,MlCrlbToff.^.5,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,BysCovToff.^.5,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,BysCrlbToff.^.5,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Error STD (Sec)'; grid on
legend('ML RMSE','ML CRLB^{1/2}',' BYS RMSE','BYS CRLB^{1/2}')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'Backup and Results\ToffStdEr_SynthData.fig')
saveas(gcf,'Backup and Results\ToffStdEr_SynthData.eps','epsc')

