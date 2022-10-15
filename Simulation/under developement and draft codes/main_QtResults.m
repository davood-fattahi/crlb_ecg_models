clc
clear
% close all


%% load the case names in the directory
DA='Backup and Results';
files=dir([DA '\*TempQt*.mat']); 
NumCases=1:2;%length(files);
i=0;
h= waitbar(i/length(NumCases), ['Case No. 0, please wait ...']);
for ii=NumCases
    i=i+1;
    waitbar(i/length(NumCases), h, ['Case No. ' num2str(ii) ', please wait ...']);

    Vars = {'SNRdB','PrParamsQ', 'PrParamsT', 'MlFINumT', 'MlFINumQ' , ...
        'MlFIApprT', 'MlFIApprQ', 'MlErQ', 'BysErQ' , 'MlErT', 'BysErT'};
    load([DA '\' files(ii).name], Vars{:}); % address of the results   
    
    if i==1
        PrParamsQ_All=PrParamsQ;
        PrParamsT_All=PrParamsT;
        MlFINumT_All=MlFINumT;
        MlFINumQ_All=MlFINumQ;
        MlFIApprT_All=MlFIApprT;
        MlFIApprQ_All=MlFIApprQ;
        MlErQ_All=MlErQ;
        MlErT_All=MlErT;
        BysErQ_All=BysErQ;
        BysErT_All=BysErT;
    else
        PrParamsQ_All=cat(1,PrParamsQ_All, PrParamsQ);
        PrParamsT_All=cat(1,PrParamsT_All, PrParamsT);
        MlFINumT_All=cat(3,MlFINumT_All,MlFINumT);
        MlFINumQ_All=cat(3,MlFINumQ_All,MlFINumQ);
        MlFIApprT_All=cat(3,MlFIApprT_All,MlFIApprT);
        MlFIApprQ_All=cat(3,MlFIApprQ_All,MlFIApprQ);
        MlErQ_All=cat(2,MlErQ_All,MlErQ);
        MlErT_All=cat(2,MlErT_All,MlErT);
        BysErQ_All=cat(2,BysErQ_All,BysErQ);
        BysErT_All=cat(2,BysErT_All,BysErT);
    end
    
end
close(h)
PrCovT_All=cov(PrParamsT_All);
PrCovQ_All=cov(PrParamsQ_All);



%%% pre-allocation
BysCrlbApprQ=zeros(3, 3, length(SNRdB));
BysCrlbApprT=zeros(3, 3, length(SNRdB));
MlCrlbApprQ_AvrgFI=zeros(3, 3, length(SNRdB));
MlCrlbApprT_AvrgdFI=zeros(3, 3, length(SNRdB));
DetMlCrlbApprQ=zeros(length(SNRdB),1);
DetBysCrlbApprQ=zeros(length(SNRdB),1);
DetMlCrlbApprT=zeros(length(SNRdB),1);
DetBysCrlbApprT=zeros(length(SNRdB),1);

BysCrlbNumQ=zeros(3, 3, length(SNRdB));
BysCrlbNumT=zeros(3, 3, length(SNRdB));
MlCrlbNumQ_AvrgFI=zeros(3, 3, length(SNRdB));
MlCrlbNumT_AvrgdFI=zeros(3, 3, length(SNRdB));
DetMlCrlbNumQ=zeros(length(SNRdB),1);
DetBysCrlbNumQ=zeros(length(SNRdB),1);
DetMlCrlbNumT=zeros(length(SNRdB),1);
DetBysCrlbNumT=zeros(length(SNRdB),1);


MlCovQ=zeros(3,3,length(SNRdB));
BysCovQ=zeros(3,3,length(SNRdB));

DetMlCovQ=zeros(size(SNRdB)); 
DetBysCovQ=zeros(size(SNRdB)); 

MlCovT=zeros(3,3,length(SNRdB));
BysCovT=zeros(3,3,length(SNRdB));

DetMlCovT=zeros(size(SNRdB)); 
DetBysCovT=zeros(size(SNRdB)); 

for k=1:length(SNRdB)

    MlCovQ(:,:,k)=cov(reshape(MlErQ_All(:,:,:,k),3,[])'); % ML Error covariance matrix
    DetMlCovQ(k)=det(MlCovQ(:,:,k)); % determinant of ML Error covariance matrix
    BysCovQ(:,:,k)=cov(reshape(BysErQ_All(:,:,:,k),3,[])'); % Bys Error covariance matrix
    DetBysCovQ(k)=det(BysCovQ(:,:,k)); % determinant of Bys Error covariance matrix

    MlCovT(:,:,k)=cov(reshape(MlErT_All(:,:,:,k),3,[])'); % ML Error covariance matrix
    DetMlCovT(k)=det(MlCovT(:,:,k)); % determinant of ML Error covariance matrix
    BysCovT(:,:,k)=cov(reshape(BysErT_All(:,:,:,k),3,[])'); % Bys Error covariance matrix
    DetBysCovT(k)=det(BysCovT(:,:,k)); % determinant of Bys Error covariance matrix

    MlCrlbNumQ_AvrgFI(:,:,k)=inv(mean(MlFINumQ_All(:,:,:,k),3,'omitnan'));
    BysCrlbNumQ(:,:,k)=inv(mean(MlFINumQ_All(:,:,:,k),3,'omitnan')+inv(PrCovQ_All));

    MlCrlbNumT_AvrgdFI(:,:,k)=inv(mean(MlFINumT_All(:,:,:,k),3,'omitnan'));
    BysCrlbNumT(:,:,k)=inv(mean(MlFINumT_All(:,:,:,k),3,'omitnan')+inv(PrCovT_All));

    DetMlCrlbNumQ(k)=det(MlCrlbNumQ_AvrgFI(:,:,k));
    DetBysCrlbNumQ(k)=det(BysCrlbNumQ(:,:,k));

    DetMlCrlbNumT(k)=det(MlCrlbNumT_AvrgdFI(:,:,k));
    DetBysCrlbNumT(k)=det(BysCrlbNumT(:,:,k));

    MlCrlbApprQ_AvrgFI(:,:,k)=inv(mean(MlFIApprQ_All(:,:,:,k),3));
    BysCrlbApprQ(:,:,k)=inv(mean(MlFIApprQ_All(:,:,:,k),3)+inv(PrCovQ_All));

    MlCrlbApprT_AvrgdFI(:,:,k)=inv(mean(MlFIApprT_All(:,:,:,k),3));
    BysCrlbApprT(:,:,k)=inv(mean(MlFIApprT_All(:,:,:,k),3)+inv(PrCovT_All));

    DetMlCrlbApprQ(k)=det(MlCrlbApprQ_AvrgFI(:,:,k));
    DetBysCrlbApprQ(k)=det(BysCrlbApprQ(:,:,k));

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
