clc
clear
close all


%% load the case names in the directory
DA='Backup and Results';
files=dir([DA '\*TempQt*.mat']); 
NumCases=length(files);

DetMlCrlbApprQ_All=0;
DetBysCrlbApprQ_All=0;
DetMlCrlbApprT_All=0;
DetBysCrlbApprT_All=0;

DetMlCrlbNumQ_All=0;
DetBysCrlbNumQ_All= 0;
DetMlCrlbNumT_All=0;
DetBysCrlbNumT_All=0;

DetMlCovQ_All=0; 
DetBysCovQ_All=0; 
DetMlCovT_All=0; 
DetBysCovT_All=0; 
for i=1:NumCases
    Vars = {'SNRdB', 'DetMlCrlbApprQ', 'DetBysCrlbApprQ', 'DetMlCrlbApprT', ...
        'DetBysCrlbApprT', 'DetMlCrlbNumQ', 'DetBysCrlbNumQ',  ...
        'DetMlCrlbNumT', 'DetBysCrlbNumT', 'DetMlCovQ', 'DetBysCovQ', 'DetMlCovT', 'DetBysCovT'};
    load([DA '\' files(i).name], Vars{:}); % address of the dataset   

    DetMlCrlbApprQ_All=DetMlCrlbApprQ_All+DetMlCrlbApprQ;
    DetBysCrlbApprQ_All=DetBysCrlbApprQ+DetBysCrlbApprQ_All;
    DetMlCrlbApprT_All=DetMlCrlbApprT+DetMlCrlbApprT_All;
    DetBysCrlbApprT_All=DetBysCrlbApprT+DetBysCrlbApprT_All;

    DetMlCrlbNumQ_All=DetMlCrlbNumQ_All+DetMlCrlbNumQ;
    DetBysCrlbNumQ_All= DetBysCrlbNumQ + DetBysCrlbNumQ_All;
    DetMlCrlbNumT_All=DetMlCrlbNumT+DetMlCrlbNumT_All;
    DetBysCrlbNumT_All=DetBysCrlbNumT+DetBysCrlbNumT_All;

    DetMlCovQ_All=DetMlCovQ+DetMlCovQ_All; 
    DetBysCovQ_All=DetBysCovQ_All+DetBysCovQ; 
    DetMlCovT_All=DetMlCovT+DetMlCovT_All; 
    DetBysCovT_All=DetBysCovT+DetBysCovT_All; 
end
DetMlCrlbApprQ_All=DetMlCrlbApprQ_All./NumCases;
DetBysCrlbApprQ_All=DetBysCrlbApprQ_All./NumCases;
DetMlCrlbApprT=DetMlCrlbApprT_All./NumCases;
DetBysCrlbApprT=DetBysCrlbApprT_All./NumCases;

DetMlCrlbNumQ_All=DetMlCrlbNumQ_All./NumCases;
DetBysCrlbNumQ_All=DetBysCrlbNumQ_All./NumCases;
DetMlCrlbNumT_All=DetMlCrlbNumT_All./NumCases;
DetBysCrlbNumT_All=DetBysCrlbNumT_All./NumCases;

DetMlCovQ_All=DetMlCovQ_All./NumCases; 
DetBysCovQ_All=DetBysCovQ_All./NumCases; 
DetMlCovT_All=DetMlCovT_All./NumCases; 
DetBysCovT_All=DetBysCovT_All./NumCases; 





    
%         figure(1)
%         semilogy(SNRdB,DetMlCrlbNumQ)
%         hold on
%         semilogy(SNRdB,DetMlCrlbApprQ)
%         xlabel 'SNR (dB)'; ylabel 'CRLB Det. (mV^2.Sec^4)';
%         legend('Num. CRLB','Appr. CRLB')
%         title 'Det. of Numeric CRLB vs. Approximated CRLB for Q wave parameters'
%         saveas(gcf,'NumApprCrlbDetQ.fig')
%         saveas(gcf,'NumApprCrlbDetQ.eps','epsc')
%         
%         
%         figure(2)
%         semilogy(SNRdB,DetMlCrlbNumT)
%         hold on
%         semilogy(SNRdB,DetMlCrlbApprT)
%         xlabel 'SNR (dB)'; ylabel 'CRLB Det. (mV^2.Sec^4)';
%         legend('Num. CRLB','Appr. CRLB')
%         title 'Det. of Numeric CRLB vs. Approximated CRLB for T wave parameters'
%         saveas(gcf,'NumApprCrlbDetT.fig')
%         saveas(gcf,'NumApprCrlbDetT.eps','epsc')
%         
        
        
        figure;
        semilogy(SNRdB,DetMlCovQ,'r-'); hold on;
        semilogy(SNRdB,DetMlCrlbNumQ,'b-')
        semilogy(SNRdB,DetBysCovQ,'m-');
        semilogy(SNRdB,DetBysCrlbNumQ,'g-');
        xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
        legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
        % title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
        saveas(gcf,'DetQ.fig')
        saveas(gcf,'DetQ.eps','epsc')
        
        
        figure;
        semilogy(SNRdB,DetMlCovT,'r-'); hold on;
        semilogy(SNRdB,DetMlCrlbNumT,'b-')
        semilogy(SNRdB,DetBysCovT,'m-');
        semilogy(SNRdB,DetBysCrlbNumT,'g-');
        xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
        legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
        % title 'Determinants of Error Cov. and CRLB matrices for T wave Gaus. parameters'
        saveas(gcf,'DetT.fig')
        saveas(gcf,'DetT.eps','epsc')

        % 

%% Ploting ...
% 
% figure
% plot(SNRdB,(BysMseQ(1,:)'-squeeze(BysCrlbQ(1,1,:)))); hold on
% plot(SNRdB,(MlMseQ(1,:)'-squeeze(MlCrlbQ_AvrgdFI(1,1,:)))); legend('BYS','ML'); title(' MSE - CRLB of Gauss. Amp.  for Q wave ');
% 
% figure
% plot(SNRdB,(BysMseQ(2,:)'-squeeze(BysCrlbQ(2,2,:)))); hold on
% plot(SNRdB,(MlMseQ(2,:)'-squeeze(MlCrlbQ_AvrgdFI(2,2,:)))); legend('BYS','ML'); title('MSE - CRLB of Gauss. Width for Q wave ');
% 
% figure
% plot(SNRdB,(BysMseQ(3,:)'-squeeze(BysCrlbQ(3,3,:)))); hold on
% plot(SNRdB,(MlMseQ(3,:)'-squeeze(MlCrlbQ_AvrgdFI(3,3,:)))); legend('BYS','ML'); title('MSE - CRLB of Gauss. Center  for Q wave ');
% 
% 
% 
% figure
% plot(SNRdB,(BysMseT(1,:)'-squeeze(BysCrlbT(1,1,:)))); hold on
% plot(SNRdB,(MlMseT(1,:)'-squeeze(MlCrlbT_AvrgdFI(1,1,:)))); legend('BYS','ML'); title('MSE - CRLB of Gauss. Amp.  for T wave ');
% 
% figure
% plot(SNRdB,(BysMseT(2,:)'-squeeze(BysCrlbT(2,2,:)))); hold on
% plot(SNRdB,(MlMseT(2,:)'-squeeze(MlCrlbT_AvrgdFI(2,2,:)))); legend('BYS','ML'); title('MSE - CRLB of Gauss. Width for T wave ');
% 
% figure
% plot(SNRdB,(BysMseT(3,:)'-squeeze(BysCrlbT(3,3,:)))); hold on
% plot(SNRdB,(MlMseT(3,:)'-squeeze(MlCrlbT_AvrgdFI(3,3,:)))); legend('BYS','ML'); title('MSE - CRLB of Gauss. Center  for T wave ');


% 
% figure
% plot(SNRdB,(DetBysCovQ(:)-DetBysCrlbQ(:))); hold on;
% plot(SNRdB,(DetMlCovQ(:)-DetMlCrlbQ(:))); legend('BYS','ML'); title('DetMSE - DetCRLB of Q wave Gauss. parameters');
% 
% figure
% plot(SNRdB,(DetBysCovT(:)-DetBysCrlbT(:))); hold on;
% plot(SNRdB,(DetMlCovT(:)-DetMlCrlbT(:))); legend('BYS','ML'); title('DetMSE - DetCRLB of T wave Gauss. parameters');

% 
% figure
% semilogy(SNRdB,DetMlCovQ); hold on;
% semilogy(SNRdB,DetMlCrlbNumQ)
% semilogy(SNRdB,DetBysCovQ);
% semilogy(SNRdB,DetBysCrlbNumQ);
% xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% % title 'Determinants of Error Cov. and CRLB matrices for Q wave Gaus. parameters'
% saveas(gcf,'DetQ.fig')
% saveas(gcf,'DetQ.eps','epsc')
% 
% 
% figure
% semilogy(SNRdB,DetMlCovT); hold on;
% semilogy(SNRdB,DetMlCrlbNumT)
% semilogy(SNRdB,DetBysCovT);
% semilogy(SNRdB,DetBysCrlbNumT);
% xlabel 'SNR (dB)'; ylabel 'Determinant (mV.Sec^2)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% % title 'Determinants of Error Cov. and CRLB matrices for T wave Gaus. parameters'
% saveas(gcf,'DetT.fig')
% saveas(gcf,'DetT.eps','epsc')
% 
% figure
% semilogy(SNRdB,MlMseQ(1,:)); hold on;
% semilogy(SNRdB,squeeze(MlCrlbQ_AvrgdFI(1,1,:)))
% semilogy(SNRdB,BysMseQ(1,:));
% semilogy(SNRdB,squeeze(BysCrlbQ(1,1,:)))
% xlabel 'SNR (dB)'; ylabel 'Square Error ((mV)^2)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'MSE vs CRLB of Gauss. Amp. for Q wave'
% saveas(gcf,'GaussAmpQ.fig')
% saveas(gcf,'GaussAmpQ.eps','epsc')
% 
% 
% figure
% semilogy(SNRdB,MlMseQ(2,:)); hold on;
% semilogy(SNRdB,squeeze(MlCrlbQ_AvrgdFI(2,2,:)))
% semilogy(SNRdB,BysMseQ(2,:));
% semilogy(SNRdB,squeeze(BysCrlbQ(2,2,:)))
% xlabel 'SNR (dB)'; ylabel 'Square Error (Sec^2)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'MSE vs CRLB of Gauss. width for Q wave'
% saveas(gcf,'GaussWidthQ.fig')
% saveas(gcf,'GaussWidthQ.eps','epsc')
% 
% 
% figure
% semilogy(SNRdB,MlMseQ(3,:)); hold on;
% semilogy(SNRdB,squeeze(MlCrlbQ_AvrgdFI(3,3,:)))
% semilogy(SNRdB,BysMseQ(3,:));
% semilogy(SNRdB,squeeze(BysCrlbQ(3,3,:)))    
% xlabel 'SNR (dB)'; ylabel 'Square Error (Sec^2)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'MSE vs CRLB of Gauss. Center for Q wave'
% saveas(gcf,'GaussCentQ.fig')
% saveas(gcf,'GaussCentQ.eps','epsc')
%   
% 
% 
% 
% figure
% semilogy(SNRdB,MlMseT(1,:)); hold on;
% semilogy(SNRdB,squeeze(MlCrlbT_AvrgdFI(1,1,:)))
% semilogy(SNRdB,BysMseT(1,:));
% semilogy(SNRdB,squeeze(BysCrlbT(1,1,:)))
% xlabel 'SNR (dB)'; ylabel 'Square Error ((mV)^2)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'MSE vs CRLB of Gauss. Amp. for T wave'
% saveas(gcf,'GaussAmpT.fig')
% saveas(gcf,'GaussAmpT.eps','epsc')
% 
% 
% figure
% semilogy(SNRdB,MlMseT(2,:)); hold on;
% semilogy(SNRdB,squeeze(MlCrlbT_AvrgdFI(2,2,:)))
% semilogy(SNRdB,BysMseT(2,:));
% semilogy(SNRdB,squeeze(BysCrlbT(2,2,:)))
% xlabel 'SNR (dB)'; ylabel 'Square Error (Sec^2)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'MSE vs CRLB of Gauss. width for T wave'
% saveas(gcf,'GaussWidthT.fig')
% saveas(gcf,'GaussWidthT.eps','epsc')
% 
% 
% figure
% semilogy(SNRdB,MlMseT(3,:)); hold on;
% semilogy(SNRdB,squeeze(MlCrlbT_AvrgdFI(3,3,:)))
% semilogy(SNRdB,BysMseT(3,:));
% semilogy(SNRdB,squeeze(BysCrlbT(3,3,:)))    
% xlabel 'SNR (dB)'; ylabel 'Square Error (Sec^2)';
% legend('ML MSE','ML CRLB',' BYS MSE','BYS CRLB')
% title 'MSE vs CRLB of Gauss. Center for T wave'
% saveas(gcf,'GaussCentT.fig')
% saveas(gcf,'GaussCentT.eps','epsc')
%   
