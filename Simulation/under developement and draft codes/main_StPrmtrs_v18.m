clc
clear
close all


%% the constants
N=200000; % the end sample of the signal loaded by rdsamp
N0=1; % the start sample of the signal loaded by rdsamp
Anttr='16a'; 

Chnl = 1; % ECG channel 
HrTr= 70/60; % Hz, threshold for heart rate, used in r peak detector
L1=60; % ms, short st length
SNRdB=(-20:20:40); % vector of SNRs in decibel
SNR=10.^(SNRdB/10); % SNRs
PolyOrder=1; % order of polynomials utilized in modeling ecg segments
w1=.75; % window length of the first median filter used in base line removing
w2=.9; % window length of the second median filter used in base line removing
PostNumBeats=[]; % max number of the beats envolved in noisy signals' parameter estimation stage
NumRuns=1; % number of repeats in noisy signals' parameter estimation stage
PriNumBeats=10000; % max number of the beats for each case envolved in prior parameter estimation stage (clean signals)
MaxNumBeatsPerCase=500;



%% load the case names in the directory
DA='..\mcode\database\long-term-st-database-1.0.0';
files=dir([DA '\*.dat']); 
NumCases=length(files);
DA='database\long-term-st-database-1.0.0';

ii=1;
AllBeats=[]; SegsST=[]; Rp=[]; JP=[]; Q=[]; T=[]; CaseNum=[]; Nrml=logical.empty;

%% stage 1- gathering the approporiat beats
for i=1:NumCases
    h = waitbar(0,['Case No. ' num2str(i) ', gathering the appropriate beats, please wait ...']);

    %% load the record
    oldFolder=cd('..\mcode'); % go to the wfdb toolbox root
    AdrsNm=[DA  '\' files(i).name(1:end-4)]; % address of the dataset
    
    % Import the annotation
    try
    [RPeaks,~,~,~,~,comments16a]=rdann(AdrsNm,Anttr,Chnl,N,N0);
    [signal,Fs,tm]=rdsamp(AdrsNm,Chnl,N,N0); % loading the data
    catch ME
        warning(['Case No. ' num2str(i) ', ' files(i).name(1:end-4) ' is skiped. See the error bellow:' newline ME.message]);
        cd(oldFolder)
        close(h);
        continue
    end
    cd(oldFolder) % go back to the main folder

    RPeaks(end)=[]; comments16a(end,:)=[];
    if MaxNumBeatsPerCase > length(RPeaks) || isempty(MaxNumBeatsPerCase)
        MaxNumBeatsPerCase=length(RPeaks);
    end
    rng(i);
    I=randperm(length(RPeaks),MaxNumBeatsPerCase);
    RPeaks=RPeaks(I);
    comments16a=comments16a(I,:);

    
    comments16a=split(comments16a,',');
    jp=RPeaks+floor(str2double(string(comments16a(:,end-3)))*Fs./1000);
    NumBeats=length(jp)-1;
    StLngth=ceil(L1*Fs/1000);

    %% baseline wandering removal,
    signal=signal(:,Chnl);    
    signal=signal-(BaseLine1(BaseLine1(signal', round(w1*Fs), 'md'), round(w2*Fs), 'mn'))';
    
    %%    
    AllBeats(end+1:end+NumBeats,1:2*floor(.5*Fs)+1)=nan;
    SegsST(end+1:end+NumBeats,1:StLngth)=nan;
    CaseNum(end+1:end+NumBeats)=i;
    Rp(end+1:end+NumBeats,1)=nan;
    JP(end+1:end+NumBeats,1)=nan;
    Nrml(end+1:end+NumBeats,1)=false;
    for j=1:NumBeats  % for each beat ...
        waitbar(j/NumBeats) % wait bar ...
        % get the beats, a half second around each r peak:
        AllBeats(ii,:)=signal(RPeaks(j)-floor(.5*Fs):RPeaks(j)+floor(.5*Fs));
        SegsST(ii,:)=signal(jp(j)+1:jp(j)+StLngth);

        Nrml(ii)=true;
        Rp(ii)=floor(.5*Fs)+1; % r peak index in each beat
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
save('Backup and Results\RealDataST_InitialPool16.mat');

%% stage 2- computing the prior

% Pre-allocating ...
PrParamsST=zeros(NumAllBeats,PolyOrder+1);

% wait bar
h = waitbar(0,'Estimating parameters for clean signals, please wait ...');
for j=1:PriNumBeats  % for each beat ...
    waitbar(j/PriNumBeats)
    % extract approx. SOI

    tm=(1:StLngth)/Fs; % time stamp; the J point is the reference time (zero).        
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
VarSigST=nan(PriNumBeats,1);
% for i=unique(CaseNum)
%     VarSigST(CaseNum==i)=var(SegsST(CaseNum==i,:),[],'all','omitnan');
% end
VarSigST(:)=var(SegsST,[],'all','omitnan');
save('Backup and Results\RealDataST_PriorSTPrmtrs18.mat');

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
VarNoiseAdd=VarSigST./SNR;

%%% 
for k=1:length(SNR) % for each SNR
    for m=1:NumRuns % for each run
        jj=0; % initializing
        for j=JJ
            waitbar(((k-1)*(m)*length(JJ)+(m-1)*length(JJ)+jj)/(length(SNR)*NumRuns*length(JJ))) % wait bar ...
            jj=jj+1; % number of beats counter

            rng(k+i+j); % fix the seed for noise generating
            SegNST=SegsST(j,:)+randn(size(SegsST(j,:))).*sqrt(VarNoiseAdd(j,k)); % adding noise to the signal;
            tm=(1:StLngth)/Fs; % time stamp; J point is the reference time (zero).        
            %%% for ML, fit gaussians on each noisy Q-segs
            MlParamsST(:,jj,m,k)=PolyFit(tm,SegNST,PolyOrder,0);
            %%% for Bayesian, fit gaussians on each noisy Q-segs
            BysParamsST(:,jj,m,k)=PolyFit(tm,SegNST,PolyOrder,0,MeanPrParamsST(:),CovPrParamsST,VarNoiseAdd(j,k));
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

MLFIST=zeros(PolyOrder+1, PolyOrder+1, length(SNRdB),PriNumBeats);
BysFIST=zeros(PolyOrder+1, PolyOrder+1, length(SNRdB),PriNumBeats);
MlCrlbST=zeros(PolyOrder+1, PolyOrder+1, length(SNRdB));
BysCrlbST=zeros(PolyOrder+1, PolyOrder+1, length(SNRdB));
DetMlCrlbST=zeros(length(SNRdB));
DetBysCrlbST=zeros(length(SNRdB));
if PolyOrder==1
    for k=1:length(SNRdB)
        for j=1:PriNumBeats
            tm=(1:StLngth)/Fs; % time stamp; the J point is the reference time (zero).        
            H=tm([1 1],:).^repmat((0:PolyOrder)',1,length(tm)); % polynomial basis
            MLFIST(:,:,k,j)=(H*H')./VarNoiseAdd(j,k); % ML Fisher Information
            BysFIST(:,:,k,j)=(H*H')./VarNoiseAdd(j,k)+inv2x2(CovPrParamsST); % Bys Fisher Information
        end
        MlCrlbST(:,:,k)=inv2x2(mean(MLFIST(:,:,k,:),4));
        BysCrlbST(:,:,k)=inv2x2(mean(BysFIST(:,:,k,:),4));
        DetMlCrlbST(k)=det(MlCrlbST(:,:,k));
        DetBysCrlbST(k)=det(BysCrlbST(:,:,k));
    end
else
    for k=1:length(SNRdB)
        for j=1:PriNumBeats
            tm=(1:StLngth)/Fs; % time stamp; the J point is the reference time (zero).        
            H=tm([1 1],:).^repmat((0:PolyOrder)',1,length(tm)); % polynomial basis
            MLFIST(:,:,k,j)=(H*H')./VarNoiseAdd(j,k); % ML Fisher Information
            BysFIST(:,:,k,j)=(H*H')./VarNoiseAdd(j,k)+inv2x2(CovPrParamsST); % Bys Fisher Information
        end
        MlCrlbST(:,:,k)=inv(mean(MLFIST(:,:,k,:),4));
        BysCrlbST(:,:,k)=inv(mean(BysFIST(:,:,k,:),4));
        DetMlCrlbST(k)=det(MlCrlbST(:,:,k));
        DetBysCrlbST(k)=det(BysCrlbST(:,:,k));
    end
end

figure
semilogy(SNRdB,DetMlCovST,'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,DetMlCrlbST,'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,DetBysCovST,'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,DetBysCrlbST,'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'DetSTparams_RealData.fig')
saveas(gcf,'DetSTparams_RealData.eps','epsc')

figure
semilogy(SNRdB,squeeze(MlCovST(1,1,:)),'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,squeeze(MlCrlbST(1,1,:)),'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,squeeze(BysCovST(1,1,:)),'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,squeeze(BysCrlbST(1,1,:)),'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'STlevel_RealData.fig')
saveas(gcf,'STlevel_RealData.eps','epsc')


figure
semilogy(SNRdB,squeeze(MlCovST(2,2,:)),'-s','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330]); hold on;
semilogy(SNRdB,squeeze(MlCrlbST(2,2,:)),'--+','LineWidth',1.5,'Color',[0 0.4470 0.7410])
semilogy(SNRdB,squeeze(BysCovST(2,2,:)),'-d','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
semilogy(SNRdB,squeeze(BysCrlbST(2,2,:)),'--x','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
saveas(gcf,'STslope_RealData.fig')
saveas(gcf,'STslope_RealData.eps','epsc')


save('Backup and Results\RealDataST_Results18.mat');



