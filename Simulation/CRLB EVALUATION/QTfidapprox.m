function [Qafp, Tafp] =QTfidapprox(ecg, Rpeak, fs)

% This function appriximately finds peaks and borders of the Q and T waves 
% in the input ecg.   


% Inputs:
% ecg: single channel ecg signal,
% R: R peak positions,
% fs: sampling frequency,

% Outputs:
% Each of the outputs is a matrix, in which the first, second, and third 
% columns respectively contain onset, peak and offset of the waves. 
% Qafp: approximated fiducal points of Q wave.
% Tafp: approximated fiducal points of T wave.


% Davood Fattahi, 5/5/2021
% fattahi.d@gmail.com




ecg=ecg(:); Rpeak=Rpeak(:);
% RR=Rpeak(2:end)-Rpeak(1:end-1);
Qafp=nan(length(Rpeak),3); Tafp=nan(length(Rpeak),3);
for i=1:length(Rpeak)
    Qafp(i,1)=Rpeak(i)-floor(.08*fs);    
    if Qafp(i,1)<=0
        Qafp(i,:)=nan;
    else
        Qafp(i,3)=Rpeak(i)-floor(.02*fs); 
        [~,I]=min(ecg(Qafp(i,1):Qafp(i,3))); Qafp(i,2)=I+Qafp(i,1)-1;
    end


    Tafp(i,3)=Rpeak(i)+floor(.5*fs);
    if Tafp(i,3)>length(ecg)
        Tafp(i,:)=nan;
    else
        Tafp(i,1)=Rpeak(i)+floor(.1*fs); 
        [~,I]=max(ecg(Tafp(i,1):Tafp(i,3))); 
        Tafp(i,2)=I+Tafp(i,1)-1;
    end
end

end

