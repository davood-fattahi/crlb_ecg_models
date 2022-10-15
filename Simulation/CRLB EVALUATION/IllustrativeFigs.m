%%
close all
clear 
load 'DataForIllustrativeFigs.mat'

%% 
%%% initial values ...
tt=(1:length(EcgBeat))/fs;

R=502;

Jp=550;
StLngth=120;

Q1=R-floor(.06*fs);
Q2=R-floor(.03*fs);
T1=R+floor(.28*fs);
T2=R+floor(.45*fs);

%%
%%% Adding noise
snrdb=15;
snr=10.^(snrdb/10);
% EcgBeatNoisy=EcgBeat+sqrt(var(EcgBeat)/snr)*randn(size(EcgBeat));
EcgBeatNoisy=EcgBeat; %+sqrt(var(EcgBeat)/snr)*randn(size(EcgBeat));

%% ST polynomials parameters 
tm=(1:StLngth)/fs;
PrParamST=PolyFit(tm,EcgBeatNoisy(Jp+1:Jp+StLngth),1,0);


figure
plot(tt,EcgBeatNoisy,'LineWidth',1,'Color',[0.3010 0.7450 0.9330])
hold on
plot(tt(Jp+1:Jp+StLngth),polyval(flip(PrParamST),tm),'LineWidth',3,'Color',[0 0 0.7410]);
xlabel 'Time (sec)'; ylabel 'Amplitude(mv)';
legend('ECG','Polynomials value');
set(gca, 'FontSize', 14)
saveas(gcf,'Backup and Results\StParamsExampleFigs.fig')
saveas(gcf,'Backup and Results\StParamsExampleFigs.eps','epsc')

%% QT parameters
options = struct('SpecifyObjectiveGradient',true);
PrParamsQ=GausFit(tt(Q1:Q2),EcgBeatNoisy(Q1:Q2),[min(EcgBeatNoisy(Q1:Q2)) (tt(Q2)-tt(Q1))/5 (tt(Q2)+tt(Q1))/2],[-inf 0 tt(Q1)],[inf inf tt(Q2)],options);
PrParamsT=GausFit(tt(T1:T2),EcgBeatNoisy(T1:T2),[max(EcgBeatNoisy(T1:T2)) (tt(T2)-tt(T1))/2 tt(Q1)],[-inf 0 tt(T1)],[inf inf tt(T2)],options);

figure;
plot(tt,EcgBeatNoisy,'LineWidth',1,'Color',[0.3010 0.7450 0.9330])
hold on
plot(tt(Q1:Q2),GausVal(tt(Q1:Q2),PrParamsQ),'LineWidth',3,'Color',[0 0 0.7410])    
plot(tt(T1:T2),GausVal(tt(T1:T2),PrParamsT),'LineWidth',3,'Color',[0 0 0.7410])    
xlabel 'Time (sec)'; ylabel 'Amplitude(mv)';
legend('ECG','Gaussians');
set(gca, 'FontSize', 14)

% saveas(gcf,'Backup and Results\QtParamsExampleFigs.fig')
saveas(gcf,'Backup and Results\QtParamsExampleFigs.eps','eps')






