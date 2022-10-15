function [P, Q, R, S, T] = ecgWavesSoI(ecg, Rpeak, RRf)

% This function appriximately finds peaks and borders of the Q and T waves 
% in the input ecg.   


% Inputs:
% ecg: single channel ecg signal,
% R: R peak positions,
% fs: sampling frequency,

% Outputs:
% Each of the outputs is a matrix, in which the first, second, and third 
% columns respectively contain onset, peak and offset of the waves. 
% Q: approximated fiducal points of Q wave.
% T: approximated fiducal points of T wave.


% Davood Fattahi, 5/5/2021
% fattahi.d@gmail.com




ecg=ecg(:); Rpeak=Rpeak(:);
RR=Rpeak(2:end)-Rpeak(1:end-1);
P=nan(length(Rpeak),3); Q=nan(length(Rpeak),3);
R=nan(length(Rpeak),3); S=nan(length(Rpeak),3);
T=nan(length(Rpeak),3); 
% Uafp=nan(length(Rpeak),3);


for i= 1: length(RR)
    P(i+1,1)=Rpeak(i+1)-floor(.4*RR(i));
    P(i+1,3)=Rpeak(i+1)-floor(.1*RR(i));
    [~,I]=max(abs(ecg(P(i+1,1):P(i+1,3)))); P(i+1,2)=I+P(i+1,1)-1;
    
    Q(i+1,1)=Rpeak(i+1)-floor(.1*RR(i));    
    Q(i+1,3)=Rpeak(i+1)-floor(.025*RR(i)); 
    [~,I]=min(ecg(Q(i+1,1):Q(i+1,3))); Q(i+1,2)=I+Q(i+1,1)-1;
    
    S(i,3)=Rpeak(i)+floor(.15*(RR(i)));
    S(i,1)=Rpeak(i)+floor(.04*(RR(i))); 
    [~,I]=min(ecg(S(i,1):S(i,3))); S(i,2)=I+S(i,1)-1;
    
    R(i,1)=Q(i,3);  R(i,3)=S(i,1); R(i,2) = Rpeak(i);

    T(i,3)=Rpeak(i)+floor(.6*RR(i));
    T(i,1)=Rpeak(i)+floor(.15*RR(i));
    [~,I]=max(ecg(T(i,1):T(i,3))); T(i,2)=I+T(i,1)-1;
end

% 
% fidp=[Q(:); R(:); S(:); T(:)]; fidp=fidp(~isnan(fidp));
% figure;
% plot(ecg); hold on; plot(fidp, ecg(fidp), 'r*');
% 


if RRf
    P=P-Rpeak; Q=Q-Rpeak; R=R-Rpeak; S=S-Rpeak; T=T-Rpeak;
end

end

