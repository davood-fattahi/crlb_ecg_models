clc
clear
close all


%% the constants
N=[]; % the end sample of the signal loaded by rdsamp
N0=1; % the start sample of the signal loaded by rdsamp
nc=20; % the n-th normal case in the database
Chnl = (1:3); % ECG channels whose mean is used as the ECG; % 2->ii, 8->v2, 11->v5
HrTr=70/60; % approximate heart rate (Hz)
NumGaus=7; % number of Gaussians utilized in modeling ecg beats
w1=.25; % window length of the first median filter used in base line removing
w2=.3; % window length of the second median filter used in base line removing
AllBeats=1000; % number of the beats envolved in signals' parameter estimation stage

%%% option structure for non linear least squre optimization, respectively
%%% in ML estimation and Bayesian estimation:
options = struct('SpecifyObjectiveGradient',true, 'FunctionTolerance', 1e-8, 'MaxFunctionEvaluations',10000, 'MaxIterations', 10000, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12);
% options = struct('SpecifyObjectiveGradient',true, 'FunctionTolerance', 1e-8, 'Display','final-detailed' ,'MaxFunctionEvaluations',10000, 'MaxIterations', 10000, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12);
% options=struct('SpecifyObjectiveGradient',true);


%% load the case names in the directory
NC=importdata('..\mcode\database\ptb-diagnostic-ecg-database-1.0.0\CONTROLS'); % load the normal cases
NumCases=length(NC);
DA='database\ptb-diagnostic-ecg-database-1.0.0';

%% load the record
oldFolder=cd('..\mcode');
AdrsNm=[DA  '\' NC{nc}];
[ecg,fs,tm]=rdsamp(AdrsNm,Chnl,N,N0); 
cd(oldFolder)

ecg=mean(ecg,2);

%% baseline wandering removal,
ecg=ecg-(BaseLine1(BaseLine1(ecg', round(w1*fs), 'md'), round(w2*fs), 'mn'))';

%% extract the R peaks using peak detector
Rp = PeakDetection(ecg(:),HrTr/fs); 
Rp=find(Rp);
%% adjusting the AllBeats by the detected R peaks:
if isempty(AllBeats) || (AllBeats >length(Rp)-2) 
    NumAllBeats=length(Rp)-2;
else
    NumAllBeats=AllBeats;
end

figure
plot(tm,ecg); hold on

% wait bar
h = waitbar(0,'Estimating parameters for clean signals, please wait ...');
prmtr=zeros(NumAllBeats,NumGaus*3);
for j=1:NumAllBeats  % for each beat ...
    waitbar(j/NumAllBeats)
    %% Q wave model fitting
    beat=ecg(floor(.4*Rp(j)+.6*Rp(j+1)):floor(.4*Rp(j+1)+.6*Rp(j+2))); % Segment Of Intrest (SOI)
    tt=tm(floor(.4*Rp(j)+.6*Rp(j+1)):floor(.4*Rp(j+1)+.6*Rp(j+2)))-tm(Rp(j)); % SOI time stamp

    % fit gssns on each beat
    prmtr(j,:)=GausFitAuto(tt,beat,NumGaus,options);

    % polt the evaluated gaussians on the signal
    ecgVal=GausVal(tt,prmtr(j,:));
    hold on;
    plot(tt+tm(Rp(j)),ecgVal,'r-')
    plot(tm(Rp(j))+prmtr(j,(1:NumGaus)*3),GausVal(prmtr(j,(1:NumGaus)*3),prmtr(j,:)),'r*');

end
close(h)



