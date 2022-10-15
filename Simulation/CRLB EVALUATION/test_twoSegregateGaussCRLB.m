clear
close all
clc



qParams=[-0.106174169797765,0.0164574890786734,-0.0475235370524113];
tParams=[0.206111783407656,0.0482330982142699,0.277807910883023];
varn = 0.1;
Ts=.001;

t=-100:Ts:100;

a1=qParams(1);
a2=tParams(1);

b1=qParams(2);
b2=tParams(2);
b=(b1+b2)/2;

c1=0:-.005:-.2;
c2=0.01:.01:0.41;
d=abs(c2-c1)/sqrt(b1.^2+b2.^2);

d_qt=[.28 .45]/sqrt(b1.^2+b2.^2);

for i=1:length(d)
    p=[a1 b1 c1(i) a2 b2 c2(i)];
    [C(:,:,i), FI(:,:,i)]=GaussCrlbNumeric(t, p, varn, 1 );
end
C = abs(C)./max(abs(C),[],3);

figure
p=plot(d,squeeze(C(1,4,:)),'LineWidth',1); hold on
plot(d,squeeze(C(2,5,:)),'LineWidth',1)
plot(d,squeeze(C(3,6,:)),'LineWidth',1)
plot(d,squeeze(C(1,5,:)),'LineWidth',1)
plot(d,squeeze(C(2,6,:)),'LineWidth',1)
plot(d,squeeze(C(1,6,:)),'LineWidth',1)
grid on

Yl=ylim;
patch([d_qt(1) d_qt(1) d_qt(2) d_qt(2)],[Yl(1) Yl(2) Yl(2) Yl(1)],'k', 'FaceAlpha', .1);
% ylim([10e-35 1])
% legend('C_{a_Qa_T} (mv^2)','C_{a_Qb_T} (mv.sec)','C_{a_Qc_T} (mv.sec)','C_{b_Qb_T} (sec^2)','C_{b_Qc_T}  (sec^2)','C_{c_Qc_T}  (sec^2)', 'd range in ECG', 'Location', 'southwest')
% legend('C_{a_1a_2}/(a_1a_2)','C_{a_1b_2}/(a_1b_2)','C_{a_1c_2}/(a_1c_2)','C_{b_1b_2}/(b_1b_2)','C_{b_1c_2}/(b_1c_2)','C_{c_1c_2}/(c_1c_2)', 'd range for Q & T', 'Location', 'southwest')
legend('\epsilon_{a_1a_2}','\epsilon_{a_1b_2}','\epsilon_{a_1c_2}','\epsilon_{b_1b_2}','\epsilon_{b_1c_2}','\epsilon_{c_1c_2}', ['d range' newline 'for Q & T'], 'Location', 'northeast')
xlabel d
ylabel 'normalized approximation error'
% ylabel('$\frac{|\epsilon_{12}|}{max(|\epsilon_{12}|)}$','fontsize', 18, 'Interpreter','latex')
xlim([0 10])

saveas(gcf,'twoGaussCrlbSimulation.fig','fig')
saveas(gcf, 'twoGaussCrlbSimulation.eps', 'epsc')

figure
p=semilogy(d,squeeze(C(1,4,:)),'LineWidth',1); hold on
semilogy(d,squeeze(C(2,5,:)),'LineWidth',1)
semilogy(d,squeeze(C(3,6,:)),'LineWidth',1)
semilogy(d,squeeze(C(1,5,:)),'LineWidth',1)
semilogy(d,squeeze(C(2,6,:)),'LineWidth',1)
semilogy(d,squeeze(C(1,6,:)),'LineWidth',1)
grid on


Yl=ylim;
patch([d_qt(1) d_qt(1) d_qt(2) d_qt(2)],[Yl(1) Yl(2) Yl(2) Yl(1)],'k', 'FaceAlpha', .1);
% legend('C_{a_1a_2}/(a_1a_2)','C_{a_1b_2}/(a_1b_2)','C_{a_1c_2}/(a_1c_2)','C_{b_1b_2}/(b_1b_2)','C_{b_1c_2}/(b_1c_2)','C_{c_1c_2}/(c_1c_2)', 'd range for Q & T', 'Location', 'southwest')
legend('\epsilon_{a_1a_2}','\epsilon_{a_1b_2}','\epsilon_{a_1c_2}','\epsilon_{b_1b_2}','\epsilon_{b_1c_2}','\epsilon_{c_1c_2}', ['d range' newline 'for Q & T'], 'Location', 'southwest')
xlabel d
ylabel 'normalized approximation error'
% ylabel('$\frac{|\epsilon_{12}|}{max(|\epsilon_{12}|)}$','fontsize', 18, 'Interpreter','latex')
xlim([0 10])

saveas(gcf,'twoGaussCrlbSimulation_log.fig','fig')
saveas(gcf, 'twoGaussCrlbSimulation_log.eps', 'epsc')


