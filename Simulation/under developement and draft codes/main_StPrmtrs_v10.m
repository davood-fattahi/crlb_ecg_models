clc
clear
close all


%% the constants
N=[];
N0=1;
Anttr='16a'; 
Chnl = 1; 
HrTr= 120; % bpm, threshold for heart rate
L1=60; % ms, short st length
L2=80; % ms, long st length
% SNR=[0.5 1 2 5 7 10:10:50];
% SNRdB=10*log10(SNR);
SNRdB=(-20:5:60); % vector
SNR=10.^(SNRdB/10);
PolyOrder=1; 
w1=.75;
w2=.9;
NumBeats=1000;
NumRuns=3;
AllBeats=[];



%% load the case names in the directory
DA='..\mcode\database\long-term-st-database-1.0.0';
files=dir([DA '\*.dat']); 
NumCases=length(files);
% NumCases=1;

DA='database\long-term-st-database-1.0.0';
for i=1:NumCases
    %% load the record
    oldFolder=cd('..\mcode');
    AdrsNm=[DA  '\' files(i).name(1:end-4)];
    [signal,Fs,~]=rdsamp(AdrsNm,[],N,N0);
    

    %% load the annotation
    % extract the J-points
    [Ann16a,~,~,~,~,comments16a]=rdann(AdrsNm,Anttr,Chnl,N,N0);
    comments16a=split(comments16a,',');
    J=Ann16a+floor(str2double(string(comments16a(:,end-3)))*Fs./1000);

    
    if isempty(AllBeats)
        NumAllBeats=length(J)-1;
    end

    cd(oldFolder)    
    %% baseline wandering removal,
    signal=signal(1:J(NumAllBeats+1),:);
    signal=signal(:,Chnl)-(BaseLine1(BaseLine1(signal(:,Chnl)', round(w1*Fs), 'md'), round(w2*Fs), 'mn'))';
 
    %% computing the prior
    % Predefine ...
    PrParams=zeros(NumAllBeats,PolyOrder+1);
    StLngth=ceil(L1*Fs/1000);
    t=(1:StLngth)*1000/Fs; tt=t-mean(t);
    StSegs=zeros(NumAllBeats,StLngth);
%     StSegEv=zeros(AllBeats,StLngth);
    % wait bar
    h = waitbar(0,['Case No. ' num2str(i) ', Estimating parameters for clean signals, please wait ...']);
    for j=1:NumAllBeats
        waitbar(j/NumAllBeats)

        %% ST segment extraction
        StSegs(j,:)=signal(J(j)+1:J(j)+StLngth);
        
        % fit plynmls on each clean ST-segs
        PrParams(j,:)=PolyFit(tt,StSegs(j,:),PolyOrder,0);
%         StSegEv(j,:)=polyval(flip(PrParams(j,:)),tt);
    end
    close(h)
    
    %%% Error histogram
%     MdlErr=StSegs-StSegEv;
%     histogram(MdlErr(:))
    
    %%% estimate the prior for params (cov matrix of params)
    PrMean=mean(PrParams,1);
    PrCov=cov(PrParams);
    VarSig=var(StSegs(:)); % variance of the signal
    
    %% add noise to the clean ST-seg
    h = waitbar(0,['Case No. ' num2str(i) ', Estimating parameters for noisy signals, please wait ...']);
    if isempty(NumBeats)
        NumBeats=NumAllBeats;
    end
    MlParams=zeros(PolyOrder+1, NumBeats, NumRuns, length(SNR));
    MlEr=zeros(size(MlParams)); 
    BysParams=zeros(PolyOrder+1, NumBeats, NumRuns, length(SNR));
    BysEr=zeros(size(BysParams));
    
    NoiseAdd=zeros(StLngth, NumBeats, NumRuns, length(SNR));
    VarNoiseAdd=zeros(size(SNR)); 
    DetMlCov=zeros(size(SNR)); 
    DetBysCov=zeros(size(SNR)); 
    
    MlCov=zeros(PolyOrder+1,PolyOrder+1,length(SNR));
    BysCov=zeros(PolyOrder+1,PolyOrder+1,length(SNR));
    
    rng(i);
    JJ=randperm(size(StSegs,1),NumBeats);
    for k=1:length(SNR)
        VarNoiseAdd(k)=VarSig/SNR(k);
        rng(k+i);
        NoiseAdd(:,:,:,k)=randn(StLngth, NumBeats, NumRuns).*sqrt(VarNoiseAdd(k));
        for m=1:NumRuns
            waitbar(((k-1)*NumRuns+m)/(length(SNR)*NumRuns))
            StSegN=StSegs(JJ,:)+NoiseAdd(:,:,m,k)';
            for j=1:NumBeats
                %%% for ML, fit plynmls on each noisy ST-segs
                MlParams(:,j,m,k)=PolyFit(tt,StSegN(j,:),PolyOrder,0);
                
                %%% for Bayesian, fit plynmls on each noisy ST-segs
                BysParams(:,j,m,k)=PolyFit(tt,StSegN(j,:),PolyOrder,0,PrMean(:),PrCov,VarNoiseAdd(k));
            end
            MlEr(:,:,m,k)=MlParams(:,:,m,k)-PrParams(JJ,:)';
            BysEr(:,:,m,k)=BysParams(:,:,m,k)-PrParams(JJ,:)';
        end
        MlCov(:,:,k)=cov(reshape(MlEr(:,:,:,k),PolyOrder+1,[])');
        DetMlCov(k)=det(MlCov(:,:,k));
        BysCov(:,:,k)=cov(reshape(BysEr(:,:,:,k),PolyOrder+1,[])');
        DetBysCov(k)=det(BysCov(:,:,k));
    end
    close(h)
    MlMse=squeeze(mean((MlParams-PrParams(JJ,:)').^2, [2 3]));
    BysMse=squeeze(mean((BysParams-PrParams(JJ,:)').^2, [2 3]));
    
    
%     save(['Backup and Results\TempStPrmtrs_' files(i).name(1:end-4)  '.mat'], 'PrParams', 'MlParams', 'BysParams', 'VarNoiseAdd', 'VarSig', 'SNR', 'MlCov', 'DetMlCov', 'BysCov', 'DetBysCov', 'MlMse', 'BysMse')    
 
clear BysEr Ann16a MlEr comments16a signal StSegs tm StSegN NoiseAdd
save(['Backup and Results\TempStPrmtrs_' files(i).name(1:end-4)  '.mat']);
end


