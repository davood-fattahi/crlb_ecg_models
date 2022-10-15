clc
clear
close all


%% the constants
Chnl = [5 2 5 11 4  4  6 2 1 2  4  4 1 7  3 1 12 1 4 1 1 1  4 6;
        1 1 1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1 1 1 1 -1 1 1 1 -1 1]; % ECG channels and their polarity.
SynthEcgFs=1000;
pstd=.03; % parameters relative standard deviation around their mean.
rrstd=.04; % rr interval relative standard deviation around its mean.
SynthSigLength=1*2*60; % length of synthetic ecg (sec)
HrTr= 70/60; % Hz, threshold for heart rate, used in r peak detector
L1=60; % ms, short st length
SNRdB=(-20:5:40); % vector of SNRs in decibel
SNR=10.^(SNRdB/10); % SNRs
PolyOrder=1; % order of polynomials utilized in modeling ecg segments
% w1=.75; % window length of the first median filter used in base line removing
% w2=.9; % window length of the second median filter used in base line removing
PostNumBeats=[]; % max number of the beats envolved in noisy signals' parameter estimation stage
NumRuns=3; % number of repeats in noisy signals' parameter estimation stage
PriNumBeats=100000; % max number of the beats for each case envolved in prior parameter estimation stage (clean signals)

%% load the case names in the directory
DA='..\mcode\database\ptb-gaussian-parameters';
files=dir([DA '\*.mat']); 
NumCases=length(files);


ii=1;
AllBeats=[]; SegsST=[]; Rp=[]; JP=[]; Q=[]; T=[]; CaseNum=[]; Nrml=logical.empty;

%% stage 1- gathering the approporiat beats
for i=1:NumCases
    h = waitbar(0,['Case No. ' num2str(i) ', gathering the appropriate beats, please wait ...']);

    %% load the parameters
    GausParam=load([DA '\' files(i).name]);
    pmean=[1e-3*Chnl(2,i).*GausParam.params{Chnl(1,i)}.a; GausParam.params{Chnl(1,i)}.b; GausParam.params{Chnl(1,i)}.theta];
    
    %% resemble the signal
    tm=1/SynthEcgFs:1/SynthEcgFs:SynthSigLength;
    rng(i) % fixing the seed for random variation in beats
    [signal, ~]=ECGResembler(tm, pmean, pstd, 1/HrTr, rrstd);
    signal=signal(:);
    StLngth=ceil(L1*SynthEcgFs/1000);

    %% extract the R peaks using peak detector
    RPeaks = PeakDetection(signal(:),HrTr/SynthEcgFs); 
    RPeaks = find(RPeaks); RPeaks([1 end])=[];
    NumBeats=length(RPeaks);
    jp=RPeaks+.040*SynthEcgFs;
    %% baseline wandering removal,
%     signal=signal(:,Chnl);    
%     signal=signal-(BaseLine1(BaseLine1(signal', round(w1*SynthEcgFs), 'md'), round(w2*SynthEcgFs), 'mn'))';
    
    %%    
    AllBeats(end+1:end+NumBeats,1:2*floor(.5*SynthEcgFs)+1)=nan;
    SegsST(end+1:end+NumBeats,1:StLngth)=nan;
    CaseNum(end+1:end+NumBeats)=i;
    Rp(end+1:end+NumBeats,1)=nan;
    JP(end+1:end+NumBeats,1)=nan;
    Nrml(end+1:end+NumBeats,1)=false;
    for j=1:NumBeats  % for each beat ...
        waitbar(j/NumBeats) % wait bar ...
        % get the beats, a half second around each r peak:
        AllBeats(ii,:)=signal(RPeaks(j)-floor(.5*SynthEcgFs):RPeaks(j)+floor(.5*SynthEcgFs));
        SegsST(ii,:)=signal(jp(j)+1:jp(j)+StLngth);

        Nrml(ii)=true;
        Rp(ii)=floor(.5*SynthEcgFs)+1; % r peak index in each beat
        JP(ii)=Rp(ii)+jp(j)-RPeaks(j); % j point index in each beat
        CaseNum(ii)=i;
        ii=ii+1;
        % Note: Q, Rp and T may be the same for all the beats, but we allocate
        % them separately for each beat for possibel changes in future. 
    end   
    close(h);
end
AllBeats=AllBeats(Nrml,:); SegsST=SegsST(Nrml,:); Rp=Rp(Nrml); JP=JP(Nrml); CaseNum=CaseNum(Nrml);
NumAllBeats=length(Rp);

%% adjusting the PriNumBeats by NumAllBeats:
if isempty(PriNumBeats) || (PriNumBeats > NumAllBeats)
    PriNumBeats=NumAllBeats;
end

%%% selecting beats randomly
rng(2); % fixing the seed for random permuting
JJ=randperm(NumAllBeats,PriNumBeats);
AllBeats=AllBeats(JJ,:); SegsST=SegsST(JJ,:); Rp=Rp(JJ); JP=JP(JJ); CaseNum=CaseNum(JJ);
save('Backup and Results\SynthDataST_InitialPool.mat');

%% stage 2- computing the prior

% Pre-allocating ...
PrParamsST=zeros(NumAllBeats,PolyOrder+1);

% wait bar
h = waitbar(0,'Estimating parameters for clean signals, please wait ...');
for j=1:PriNumBeats  % for each beat ...
    waitbar(j/PriNumBeats)
    % extract approx. SOI

    tm=((1:StLngth)-mean(1:StLngth))/SynthEcgFs; % time stamp; SoI's center is the reference time (zero).        
    PrParamsST(j,:)=PolyFit(tm,SegsST(j,:),PolyOrder,0);

    %         StSegEv(j,:)=polyval(flip(PrParams(j,:)),tt);

%     % polt the evaluated gaussians on the signal
%     hold off; plot(tm,SegsST(j,:)); hold on
%     plot(tm,polyval(flip(PrParams(j,:)),tm),'r-')
end
close(h)

% estimate the prior for the params
MeanPrParamsST=mean(PrParamsST,1);
CovPrParamsST=cov(PrParamsST);
% VarSigST=nan(PriNumBeats,1);
% for i=unique(CaseNum)
%     VarSigST(CaseNum==i)=var(SegsST(CaseNum==i,:),[],'all','omitnan');
% end
save('Backup and Results\SynthDataSTcenteredSoI_PriorSTPrmtrs.mat');

%% Stage 3- estimating the parameters from noisy signals
h = waitbar(0,'Estimating parameters for noisy signals, please wait ...');
%%% adjusting the number of noisy beats
if isempty(PostNumBeats)
    PostNumBeats=PriNumBeats;
end
if PostNumBeats > PriNumBeats
    PostNumBeats=PriNumBeats;
end

%%% pre allocating ...
MlParamsST=zeros(PolyOrder+1, PostNumBeats, NumRuns, length(SNR));
MlErST=zeros(size(MlParamsST));     

BysParamsST=zeros(PolyOrder+1, PostNumBeats, NumRuns, length(SNR));
BysErST=zeros(size(BysParamsST));     

MlCovST=zeros(PolyOrder+1,PolyOrder+1,length(SNR));
BysCovST=zeros(PolyOrder+1,PolyOrder+1,length(SNR));

DetMlCovST=zeros(size(SNR)); 
DetBysCovST=zeros(size(SNR)); 

 
%%% select the noisy beats randomly
rng(3) % fixing the seed for random permuting
JJ=randperm(PriNumBeats,PostNumBeats);
VarNoiseAdd=var(SegsST,[],'all')./SNR;

%%% 
for k=1:length(SNR) % for each SNR
    for m=1:NumRuns % for each run
        jj=0; % initializing
        for j=JJ
            waitbar(((k-1)*(NumRuns)*length(JJ)+(m-1)*length(JJ)+jj)/(length(SNR)*NumRuns*length(JJ))) % wait bar ...
            jj=jj+1; % number of beats counter

            rng(k+i+j); % fix the seed for noise generating
            SegNST=SegsST(j,:)+randn(size(SegsST(j,:))).*sqrt(VarNoiseAdd(k)); % adding noise to the signal;
            tm=((1:StLngth)-mean(1:StLngth))/SynthEcgFs; % time stamp; SoI's center is the reference time (zero).        
            %%% for ML, fit gaussians on each noisy Q-segs
            MlParamsST(:,jj,m,k)=PolyFit(tm,SegNST,PolyOrder,0);
            %%% for Bayesian, fit gaussians on each noisy Q-segs
            BysParamsST(:,jj,m,k)=PolyFit(tm,SegNST,PolyOrder,0,MeanPrParamsST(:),CovPrParamsST,VarNoiseAdd(k));
        end
        MlErST(:,:,m,k)=MlParamsST(:,:,m,k)-PrParamsST(JJ,:)'; % ML Error from prior
        BysErST(:,:,m,k)=BysParamsST(:,:,m,k)-PrParamsST(JJ,:)'; % Bys Error from prior
    end
    MlCovST(:,:,k)=cov(reshape(MlErST(:,:,:,k),PolyOrder+1,[])'); % ML Error covariance matrix
    DetMlCovST(k)=det(MlCovST(:,:,k)); % determinant of ML Error covariance matrix
    BysCovST(:,:,k)=cov(reshape(BysErST(:,:,:,k),PolyOrder+1,[])'); % Bys Error covariance matrix
    DetBysCovST(k)=det(BysCovST(:,:,k)); % determinant of Bys Error covariance matrix
end
close(h)

%% CRLB calculation
% polynomials bases
H=tm([1 1],:).^repmat((0:PolyOrder)',1,length(tm)); % polynomial basis

MlCrlbST=zeros(PolyOrder+1, PolyOrder+1, length(SNRdB));
BysCrlbST=zeros(PolyOrder+1, PolyOrder+1, length(SNRdB));
DetMlCrlbST=zeros(length(SNRdB),1);
DetBysCrlbST=zeros(length(SNRdB),1);
if PolyOrder==1
    for k=1:length(SNRdB)
        MlCrlbST(:,:,k)=inv2x2((H*H')./VarNoiseAdd(k));
        BysCrlbST(:,:,k)=inv2x2((H*H')./VarNoiseAdd(k)+inv2x2(CovPrParamsST));
        DetMlCrlbST(k)=det(MlCrlbST(:,:,k));
        DetBysCrlbST(k)=det(BysCrlbST(:,:,k));
    end
else
    for k=1:length(SNRdB)
        MlCrlbST(:,:,k)=inv((H*H')./VarNoiseAdd(k));
        BysCrlbST(:,:,k)=inv((H*H')./VarNoiseAdd(k)+inv(CovPrParamsST));
        DetMlCrlbST(k)=det(MlCrlbST(:,:,k));
        DetBysCrlbST(k)=det(BysCrlbST(:,:,k));
    end
end

save('Backup and Results\SynthDataSTcenteredSoI_Results.mat');

%% ST index error and CRLB
clear
close all
load('Backup and Results\results heavy run\SynthDataSTCenteredSoI_Results.mat');
beta=.1;

figure
semilogy(SNRdB,DetMlCovST,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,DetMlCrlbST,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,DetBysCovST,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,DetBysCrlbST,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Determinant (mV^4.Sec^{-2})'; grid on
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'Backup and Results\DetSTparams_SynthDataCenteredSoI.fig')
saveas(gcf,'Backup and Results\DetSTparams_SynthDataCenteredSoI.eps','epsc')

figure
semilogy(SNRdB,squeeze(MlCovST(1,1,:)),'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,squeeze(MlCrlbST(1,1,:)),'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,squeeze(BysCovST(1,1,:)),'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,squeeze(BysCrlbST(1,1,:)),'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Error varinace (mV^2)'; grid on
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'Backup and Results\STlevelVarEr_SynthDataCenteredSoI.fig')
saveas(gcf,'Backup and Results\STlevelVarEr_SynthDataCenteredSoI.eps','epsc')


figure
semilogy(SNRdB,squeeze(MlCovST(2,2,:)),'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,squeeze(MlCrlbST(2,2,:)),'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,squeeze(BysCovST(2,2,:)),'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,squeeze(BysCrlbST(2,2,:)),'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Error varinace (mV^2.Sec^{-2})'; grid on
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'Backup and Results\STslopeVarEr_SynthDataCenteredSoI.fig')
saveas(gcf,'Backup and Results\STslopeVarEr_SynthDataCenteredSoI.eps','epsc')


figure
semilogy(SNRdB,squeeze(MlCovST(1,1,:)),'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,squeeze(MlCrlbST(1,1,:)),'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,squeeze(BysCovST(1,1,:)),'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,squeeze(BysCrlbST(1,1,:)),'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Error STD (mV)'; grid on
legend('ML RMSE','ML CRLB^{1/2}',' BYS RMSE','BYS CRLB^{1/2}')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'Backup and Results\STlevelStdEr_SynthDataCenteredSoI.fig')
saveas(gcf,'Backup and Results\STlevelStdEr_SynthDataCenteredSoI.eps','epsc')


figure
semilogy(SNRdB,squeeze(MlCovST(2,2,:)),'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,squeeze(MlCrlbST(2,2,:)),'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,squeeze(BysCovST(2,2,:)),'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,squeeze(BysCrlbST(2,2,:)),'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Error STD (mV.Sec^{-1})'; grid on
legend('ML RMSE','ML CRLB^{1/2}',' BYS RMSE','BYS CRLB^{1/2}')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'Backup and Results\STslopeStdEr_SynthDataCenteredSoI.fig')
saveas(gcf,'Backup and Results\STslopeStdEr_SynthDataCenteredSoI.eps','epsc')

% ST index error and CRLB
PrStInd=PrParamsST(:,1)+beta*PrParamsST(:,2);
MlParamsStInd=squeeze(MlParamsST(1,:,:,:)+beta*MlParamsST(2,:,:,:));
BysParamsStInd=squeeze(BysParamsST(1,:,:,:)+beta*BysParamsST(2,:,:,:));
MlErStInd=MlParamsStInd-PrStInd(JJ); % ML Error from prior
BysErStInd=BysParamsStInd-PrStInd(JJ); % Bys Error from prior

MlCovStInd=zeros(size(SNR)); BysCovStInd=zeros(size(SNR)); MlCrlbStInd=zeros(size(SNR)); BysCrlbStInd=zeros(size(SNR));
for k=1:length(SNR)
    MlCovStInd(k)=var(MlErStInd(:,:,k),[],'all'); % ML Error variance
    BysCovStInd(k)=var(BysErStInd(:,:,k),[],'all'); % Bys Error varinace
    MlCrlbStInd(k)=[1 beta]*MlCrlbST(:,:,k)*[1 beta]';
    BysCrlbStInd(k)=[1 beta]*BysCrlbST(:,:,k)*[1 beta]';
end


figure
semilogy(SNRdB,MlCovStInd,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,MlCrlbStInd,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,BysCovStInd,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,BysCrlbStInd,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Error varinace ((mV + mV/Sec)^2)'; grid on
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'Backup and Results\STindVarEr_SynthDataCenteredSoI.fig')
saveas(gcf,'Backup and Results\STindVarEr_SynthDataCenteredSoI.eps','epsc')



figure
semilogy(SNRdB,MlCovStInd,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,MlCrlbStInd,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,BysCovStInd,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,BysCrlbStInd,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Error STD (mV+mV.Sec^{-1})'; grid on
legend('ML RMSE','ML CRLB^{1/2}',' BYS RMSE','BYS CRLB^{1/2}')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'Backup and Results\STindStdEr_SynthDataCenteredSoI.fig')
saveas(gcf,'Backup and Results\STindStdEr_SynthDataCenteredSoI.eps','epsc')

