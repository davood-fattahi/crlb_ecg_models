clear
clc
close all
% load 'signal2.mat'

%% load the case names in the directory
N=[]; % the end sample of the signal loaded by rdsamp
N0=1; % the start sample of the signal loaded by rdsamp
DA='..\mcode\database\qt-database-1.0.0';
files=dir([DA '\*.dat']); 
NumCases=length(files);
DA='database\qt-database-1.0.0';


%% load the record
oldFolder=cd('..\mcode'); % go to the wfdb toolbox root
AdrsNm=[DA  '\' files(35).name(1:end-4)]; % address of the dataset
[signal,Fs,~]=rdsamp(AdrsNm,[],N,N0); % loading the data
cd(oldFolder) % go back to the main folder
fs=Fs;


%%
ecg=signal(:,1);
w1=.6;
w2=.7;
ecgBR=ecg-(BaseLine1(BaseLine1(ecg', round(w1*fs), 'md'), round(w2*fs), 'mn'))';
ecgB=ecgBeats(ecgBR,fs);
plot((1:size(ecgB,2))/fs,ecgB(1:100,:)')
xlabel 'Time (s)'; ylabel 'Amplitude(mv)';
% legend('ECG');
set(gca, 'FontSize', 14)
grid minor

saveas(gcf,'ecgBeatsDevExample.fig')
saveas(gcf,'ecgBeatsDevExample.png')
%%

ecg=signal(1:600,1);
w1=.6;
w2=.7;
ecgBR=ecg-(BaseLine1(BaseLine1(ecg', round(w1*fs), 'md'), round(w2*fs), 'mn'))';
[ecgM, ecgMed, rp, rps]=ecgmean(ecgBR,fs);
% figure
% subplot(211)
% plot(ecgM); hold on
% plot(ecgMed);
%
% fs=Fs;
% ecg=signal(1:5000,2);
% w1=.6;
% w2=.8;
% ecgBR=ecg-(BaseLine1(BaseLine1(ecg', round(w1*fs), 'md'), round(w2*fs), 'mn'))';
% [ecgM, ecgMed, rp, rps]=ecgmean(ecgBR,fs);
% subplot(212)
% plot(ecgM); hold on
% plot(ecgMed);

%%
tt=(1:length(ecgM))/fs;
R=88;
Q1=R-floor(.06*fs);
Q2=R-floor(.03*fs);
T1=R+floor(.28*fs);
T2=R+floor(.5*fs);

%% T-off annotation 
Toffi=215;

figure
plot(tt,ecgM,'LineWidth',1,'Color',[0.3010 0.7450 0.9330])
hold on
plot(tt(Toffi),0,'ok','LineWidth',2)
aa=annotation('textarrow',[0.685 0.685],[.47 0.28],'String','T-wave offset');

aa=annotation('doublearrow', [0.53 0.68], [.22 0.22]);
aa=annotation('arrow', [ 0.83 0.69], [.22 0.22]);
aa=annotation('arrow', [ 0.13 0.28], [.22 0.22]);
text(.8,-.27,'T-wave','HorizontalAlignment','right')
text(1.1,-.27,'iso-electric','HorizontalAlignment','right')
text(.2,-.27,'iso-electric','HorizontalAlignment','right')
xlabel 'Time (s)'; ylabel 'Amplitude(mv)';
legend('ECG');

set(gca, 'FontSize', 14)

grid minor

saveas(gcf,'T0ffExampleAnnotate.fig')
saveas(gcf,'T0ffExampleAnnotate.png')

%% T-off 1-st degree polynomials parameters 
soi=T2-floor(.1*fs)+1:T2;
PParam=PolyFit(tt(soi),ecgM(soi),1,0);
Toffi=211;
figure
plot(tt,ecgM,'LineWidth',1,'Color',[0.3010 0.7450 0.9330])
hold on
plot(tt(soi),polyval(flip(PParam),tt(soi)),'LineWidth',2,'Color',[0 0 0.7410]);
plot(tt,zeros(size(tt)),'r','LineWidth',2)
plot(tt(Toffi),0,'o','LineWidth',2)
xlabel 'Time (s)'; ylabel 'Amplitude(mv)';
legend('ECG','1^{st} order polynomial', 'isoelectric');
set(gca, 'FontSize', 14)
% a=annotation('textarrow',[0.8 0.69],[.5 0.28],'String','T-wave \newline offset ');
aa=annotation('textarrow',[0.68 0.68],[.47 0.28],'String','T-wave offset, \newline intersection with isoelectric');
grid minor

saveas(gcf,'T0ffExample1poly.fig')
saveas(gcf,'T0ffExample1poly.png')



%% T-off 2-nd polynomials parameters 
soi=T2-floor(.07*fs)+1:T2+floor(.03*fs);
PParam=PolyFit(tt(soi),ecgM(soi),2,0);
Toffi=213;
figure
plot(tt,ecgM,'LineWidth',1,'Color',[0.3010 0.7450 0.9330])
hold on
plot(tt(soi),polyval(flip(PParam),tt(soi)),'LineWidth',2,'Color',[0 0 0.7410]);
plot(tt,zeros(size(tt)),'r','LineWidth',2)
plot(tt(Toffi),0,'o','LineWidth',2)

xlabel 'Time (s)'; ylabel 'Amplitude(mv)';
legend('ECG','2^{nd} order polynomial', 'isoelectric');
set(gca, 'FontSize', 14)
% a=annotation('textarrow',[0.8 0.69],[.5 0.28],'String','T-wave \newline offset ');
aa=annotation('textarrow',[0.68 0.68],[.43 0.28],'String','T-wave offset, \newline intersection with isoelectric');
grid minor
saveas(gcf,'T0ffExamp2eAnnotate.fig')
saveas(gcf,'T0ffExamp2eAnnotate.png')

%% T-off 2-nd polynomials parameters 
soi=T2-floor(.15*fs)+1:T2-floor(.01*fs);
PParam=PolyFit(tt(soi),ecgM(soi),2,0);
Toffi=209;
figure
plot(tt,ecgM,'LineWidth',1,'Color',[0.3010 0.7450 0.9330])
hold on
plot(tt(soi),polyval(flip(PParam),tt(soi)),'LineWidth',2,'Color',[0 0 0.7410]);
plot(tt,zeros(size(tt)),'r','LineWidth',2)
plot(tt(Toffi),0,'o','LineWidth',2)

xlabel 'Time (s)'; ylabel 'Amplitude(mv)';
legend('ECG','2^{nd} order polynomial', 'isoelectric');
set(gca, 'FontSize', 14)
% a=annotation('textarrow',[0.8 0.69],[.5 0.28],'String','T-wave \newline offset ');
aa=annotation('textarrow',[0.68 0.68],[.43 0.28],'String','T-wave offset, \newline intersection with isoelectric');
grid minor

% saveas(gcf,'Backup and Results\StParamsExampleFigs.fig')
% saveas(gcf,'Backup and Results\StParamsExampleFigs.eps','epsc')

%% T-off 3-rd polynomials parameters 
soi=T2-floor(.13*fs)+1:T2+floor(.01*fs);
PParam=PolyFit(tt(soi),ecgM(soi),3,0);

Toffi=212;
figure
plot(tt,ecgM,'LineWidth',1,'Color',[0.3010 0.7450 0.9330])
hold on
plot(tt(soi),polyval(flip(PParam),tt(soi)),'LineWidth',2,'Color',[0 0 0.7410]);
% plot(tt,zeros(size(tt)),'r','LineWidth',2)
plot(tt(Toffi),ecgM(Toffi),'ro','LineWidth',2)
xlabel 'Time (s)'; ylabel 'Amplitude(mv)';
legend('ECG','3^{rd} order polynomial');
set(gca, 'FontSize', 14)
% a=annotation('textarrow',[0.8 0.69],[.5 0.28],'String','T-wave \newline offset ');
aa=annotation('textarrow',[0.68 0.68],[.43 0.28],'String','T-wave offset, \newline 99% drop in amplitude ');
grid minor
saveas(gcf,'T0ffExamp3eAnnotate.fig')
saveas(gcf,'T0ffExamp3eAnnotate.png')

%% T-off Gaussian function parameters 
soi=T1:T2+floor(.06*fs);
options = struct('SpecifyObjectiveGradient',true);
PrParamsT=GausFit(tt(soi),ecgM(soi),[min(ecgM(soi)) (tt(T2)-tt(T1))/5 (tt(T2)+tt(T1))/2],[-inf 0 tt(T1)],[inf inf tt(T2)],options);
Toff=PrParamsT(3) + 3*PrParamsT(2);
Toffi=floor(Toff*fs);
figure;
plot(tt,ecgM,'LineWidth',1,'Color',[0.3010 0.7450 0.9330])
hold on
plot(tt(soi),GausVal(tt(soi),PrParamsT),'LineWidth',2,'Color',[0 0 0.7410])  
plot(tt(Toffi),ecgM(Toffi),'ro','LineWidth',2)
xlabel 'Time (sec)'; ylabel 'Amplitude(mv)';
legend('ECG','Gaussian function');
set(gca, 'FontSize', 14)
grid minor

% a=annotation('textarrow',[0.8 0.7],[.5 0.28],'String','T-wave \newline offset ');
aa=annotation('textarrow',[0.7 0.7],[.43 0.28],'String','T-wave offset, \newline 99% drop in amplitude ');
saveas(gcf,'T0ffExampleGauss.fig')
saveas(gcf,'T0ffExampleGauss.png')

% a.Color = 'red'







