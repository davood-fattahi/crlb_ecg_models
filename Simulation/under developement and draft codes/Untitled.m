clc
clear
close all


%% the constants
N=[]; % the end sample of the signal loaded by rdsamp
N0=1; % the start sample of the signal loaded by rdsamp
Chnl = 1; % ECG channel 
phase0 = 0;   % desired phase shift
bins = 250;   % number of phase bins
SynthEcgFs=1000;
pstd=.04; % parameters relative standard deviation around their mean.
rrstd=.04; % rr interval relative standard deviation around its mean.
SynthSigLength=1*1*60; % length of synthetic ecg (sec)
HrTr= 70/60; % Hz, threshold for heart rate, used in r peak detector
SNRdB=(-20:20:40); % vector of SNRs in decibel
SNR=10.^(SNRdB/10); % SNRs
NumGaus=1; % number of Gaussians utilized in modeling ecg segments
w1=.25; % window length of the first median filter used in base line removing
w2=.30; % window length of the second median filter used in base line removing
PostNumBeats=[]; % max number of the beats envolved in noisy signals' parameter estimation stage
NumRuns=1; % number of repeats in noisy signals' parameter estimation stage
PriNumBeats=[]; % max number of the beats envolved in prior parameter estimation stage (clean signals)

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
DA='..\mcode\database\ptb-gaussian-parameters';
files=dir([DA '\*.mat']); 
NumCases=length(files);

ii=1;
AllBeats=[]; Rp=[]; Q=[]; T=[]; CaseNum=[]; Nrml=[];

%% stage 1- gathering the approporiat beats
for i=1:NumCases
    %% load the parameters
    GausParam=load([DA '\' files(i).name]);
    figure('Units', 'normalized', 'position', [.05 .05  .6 .8]);
    for Chnl=1:12
        pmean=[GausParam.params{Chnl}.a; GausParam.params{Chnl}.b; GausParam.params{Chnl}.theta];
        subplot(6,2,Chnl)
        plot(GausVal(-pi:2*pi/200:pi,pmean));
    end
end