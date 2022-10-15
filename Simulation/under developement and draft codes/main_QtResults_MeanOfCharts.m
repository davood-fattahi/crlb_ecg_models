clc
clear
close all


%% load the case names in the directory
DA='Backup and Results';
files=dir([DA '\*TempQt*.mat']); 
NumCases=length(files);

MlCrlbNumQ_AvrgFI_All=0;
BysCrlbNumQ_All=0;
MlCrlbNumT_AvrgdFI_All=0;
BysCrlbNumT_All=0;
MlCrlbApprQ_AvrgFI_All=0;
BysCrlbApprQ_All=0;
MlCrlbApprT_AvrgdFI_All=0;
BysCrlbApprT_All=0;
MlCovQ_All=0;
BysCovQ_All=0;
MlCovT_All=0;
BysCovT_All=0;


for i=1:NumCases
    Vars = {'SNRdB', 'MlCovQ', 'BysCovQ', 'MlCovT', 'BysCovT', 'MlCrlbNumQ_AvrgFI', 'BysCrlbNumQ', 'MlCrlbNumT_AvrgdFI', ...
        'BysCrlbNumT', 'MlCrlbApprQ_AvrgFI', 'BysCrlbApprQ', 'MlCrlbApprT_AvrgdFI', 'BysCrlbApprT'};
    load([DA '\' files(i).name], Vars{:}); % address of the dataset   
    
    MlCrlbNumQ_AvrgFI_All=MlCrlbNumQ_AvrgFI+MlCrlbNumQ_AvrgFI_All;
    BysCrlbNumQ_All=BysCrlbNumQ+BysCrlbNumQ_All;
    MlCrlbNumT_AvrgdFI_All=MlCrlbNumT_AvrgdFI+MlCrlbNumT_AvrgdFI_All;
    BysCrlbNumT_All=BysCrlbNumT+BysCrlbNumT_All;
    MlCrlbApprQ_AvrgFI_All=MlCrlbApprQ_AvrgFI+MlCrlbApprQ_AvrgFI_All;
    BysCrlbApprQ_All=BysCrlbApprQ_All+BysCrlbApprQ;
    MlCrlbApprT_AvrgdFI_All=MlCrlbApprT_AvrgdFI_All+MlCrlbApprT_AvrgdFI;
    BysCrlbApprT_All=BysCrlbApprT+BysCrlbApprT_All;
    MlCovQ_All=MlCovQ_All+MlCovQ;
    BysCovQ_All=BysCovQ_All+BysCovQ;
    MlCovT_All=MlCovT_All+MlCovT;
    BysCovT_All=BysCovT_All+BysCovT;    

end

MlCrlbNumQ_AvrgFI_All=MlCrlbNumQ_AvrgFI_All./NumCases;
BysCrlbNumQ_All=BysCrlbNumQ_All./NumCases;
MlCrlbNumT_AvrgdFI_All=MlCrlbNumT_AvrgdFI_All./NumCases;
BysCrlbNumT_All=BysCrlbNumT_All./NumCases;
MlCrlbApprQ_AvrgFI_All=MlCrlbApprQ_AvrgFI_All./NumCases;
BysCrlbApprQ_All=BysCrlbApprQ_All./NumCases;
MlCrlbApprT_AvrgdFI_All=MlCrlbApprT_AvrgdFI_All./NumCases;
BysCrlbApprT_All=BysCrlbApprT_All./NumCases;
MlCovQ_All=MlCovQ_All./NumCases;
BysCovQ_All=BysCovQ_All./NumCases;
MlCovT_All=MlCovT_All./NumCases;
BysCovT_All=BysCovT_All./NumCases;


%%% pre-allocation

DetMlCrlbApprQ=zeros(length(SNRdB),1);
DetBysCrlbApprQ=zeros(length(SNRdB),1);
DetMlCrlbApprT=zeros(length(SNRdB),1);
DetBysCrlbApprT=zeros(length(SNRdB),1);

DetMlCrlbNumQ=zeros(length(SNRdB),1);
DetBysCrlbNumQ=zeros(length(SNRdB),1);
DetMlCrlbNumT=zeros(length(SNRdB),1);
DetBysCrlbNumT=zeros(length(SNRdB),1);

DetMlCovQ=zeros(size(SNRdB)); 
DetBysCovQ=zeros(size(SNRdB)); 
DetMlCovT=zeros(size(SNRdB)); 
DetBysCovT=zeros(size(SNRdB)); 

for k=1:length(SNRdB)
    DetMlCovQ(k)=det(MlCovQ_All(:,:,k)); % determinant of ML Error covariance matrix
    DetBysCovQ(k)=det(BysCovQ_All(:,:,k)); % determinant of Bys Error covariance matrix

    DetMlCovT(k)=det(MlCovT_All(:,:,k)); % determinant of ML Error covariance matrix
    DetBysCovT(k)=det(BysCovT_All(:,:,k)); % determinant of Bys Error covariance matrix

    DetMlCrlbNumQ(k)=det(MlCrlbNumQ_AvrgFI_All(:,:,k));
    DetBysCrlbNumQ(k)=det(BysCrlbNumQ_All(:,:,k));

    DetMlCrlbNumT(k)=det(MlCrlbNumT_AvrgdFI_All(:,:,k));
    DetBysCrlbNumT(k)=det(BysCrlbNumT_All(:,:,k));

    DetMlCrlbApprQ(k)=det(MlCrlbApprQ_AvrgFI_All(:,:,k));
    DetBysCrlbApprQ(k)=det(BysCrlbApprQ_All(:,:,k));

    DetMlCrlbApprT(k)=det(MlCrlbApprT_AvrgdFI(:,:,k));
    DetBysCrlbApprT(k)=det(BysCrlbApprT(:,:,k));
end






    
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
