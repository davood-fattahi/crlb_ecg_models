
%% polynomials bases
H=zeros(PolyOrder+1,length(tt));
for i=0:PolyOrder
    H(i+1,:)=tt.^i;
end


MlCrlb=zeros(PolyOrder+1, PolyOrder+1, length(SNRdB));
BysCrlb=zeros(PolyOrder+1, PolyOrder+1, length(SNRdB));
if PolyOrder==1
    for k=1:length(SNRdB)
        MlCrlb(:,:,k)=inv2x2((H*H')./VarNoiseAdd(k));
        BysCrlb(:,:,k)=inv2x2((H*H')./VarNoiseAdd(k)+inv2x2(PrCov));
        DetMlCrlb(k)=det(MlCrlb(:,:,k));
        DetBysCrlb(k)=det(BysCrlb(:,:,k));
        EigBysCrlb(:,k)=eig(BysCov(:,:,k)-BysCrlb(:,:,k));
        EigMlCrlb(:,k)=eig(BysCov(:,:,k)-BysCrlb(:,:,k));
    end
else
    for k=1:length(SNRdB)
        MlCrlb(:,:,k)=inv((H*H')./VarNoiseAdd(k));
        BysCrlb(:,:,k)=inv((H*H')./VarNoiseAdd(k)+inv(PrCov));
        DetMlCrlb(k)=det(MlCrlb(:,:,k));
        DetBysCrlb(k)=det(BysCrlb(:,:,k));
        EigBysCrlb(:,k)=eig(BysCov(:,:,k)-BysCrlb(:,:,k));
        EigMlCrlb(:,k)=eig(MlCov(:,:,k)-MlCrlb(:,:,k));
        
    end
end
  




%% Ploting ...
% 
% figure
% plot(SNRdB,(BysMse(1,:)'-squeeze(BysCrlb(1,1,:)))); hold on
% plot(SNRdB,(MlMse(1,:)'-squeeze(MlCrlb(1,1,:)))); legend('BYS','ML'); title('MSE - CRLB of ST level');
% 
% figure
% plot(SNRdB,(BysMse(2,:)'-squeeze(BysCrlb(2,2,:)))); hold on
% plot(SNRdB,(MlMse(2,:)'-squeeze(MlCrlb(2,2,:)))); legend('BYS','ML'); title('MSE - CRLB of ST slope');
% 
% 
% figure
% plot(SNRdB,(DetBysCov-DetBysCrlb)); hold on;
% plot(SNRdB,(DetMlCov-DetMlCrlb)); legend('BYS','ML'); title('DetMSE - DetCRLB of ST polynomials coeffs');


figure
semilogy(SNRdB,DetMlCov); hold on;
semilogy(SNRdB,DetMlCrlb)
semilogy(SNRdB,DetBysCov);
semilogy(SNRdB,DetBysCrlb);
xlabel 'SNR'; ylabel 'Determinant';
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
title 'Det. of Error Cov. Matrix vs. Det. of CRLB'
saveas(gcf,'DetStLvl.fig')
saveas(gcf,'DetStLvl.eps','epsc')

figure
semilogy(SNRdB,MlMse(1,:)); hold on;
semilogy(SNRdB,squeeze(MlCrlb(1,1,:)))
semilogy(SNRdB,BysMse(1,:));
semilogy(SNRdB,squeeze(BysCrlb(1,1,:)))
xlabel 'SNR'; ylabel 'Square Error';
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
title 'ST LEVEL'
saveas(gcf,'StLvl.fig')
saveas(gcf,'StLvl.eps','epsc')


figure
semilogy(SNRdB,MlMse(2,:)); hold on;
semilogy(SNRdB,squeeze(MlCrlb(2,2,:)))
semilogy(SNRdB,BysMse(2,:));
semilogy(SNRdB,squeeze(BysCrlb(2,2,:)))    
xlabel 'SNR'; ylabel 'Square Error';
legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
title 'ST SLOPE'
saveas(gcf,'StSlp.fig')
saveas(gcf,'StSlp.eps','epsc')
    

    
    