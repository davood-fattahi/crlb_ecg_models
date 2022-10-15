function J=AnJ4LTST(AnAdrsNm, Anttr, Chnl, N, N0, Fs)

[Ann16a,~,~,~,~,comments16a]=rdann(AnAdrsNm,Anttr,Chnl,N,N0);
comments16a=split(comments16a,',');
% comments16a1=comments16a(chan16a==Chnl,:);
% Ann16a1=Ann16a(chan16a==Chnl);
J=Ann16a+floor(str2double(string(comments16a(:,end-3)))*Fs./1000);
