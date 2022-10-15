function [Pafp, Qafp, Rafp, Safp, Tafp] =ecgfidapprox2(ecg, Rpeak)

% This function appriximately finds peaks and borders of the input ecg beat
% waves.   

% [Pafp, Qafp, Rafp, Safp, Tafp] =ecgfidapprox2(ecg, Rpeak, varargin)

% Inputs:
% ecg: single channel ecg
% R: R peak positions

% Outputs:
% Each of the outputs is a matrix, in which the first, second, and third 
% columns respectively contain onset, peak and offset of the waves. 
% Pafp: approximated fiducal points of P wave.


% Davood Fattahi, 3/9/2021
% fattahi.d@gmail.com




ecg=ecg(:); Rpeak=Rpeak(:);
RR=Rpeak(2:end)-Rpeak(1:end-1);
Rafp=[nan(size(Rpeak)) Rpeak nan(size(Rpeak))];
Qafp=nan(size(Rafp)); Safp=nan(size(Rafp)); Pafp=nan(size(Rafp)); Tafp=nan(size(Rafp));
for i=2:length(Rpeak)-1
    
    Pafp(i,1)=Rpeak(i)-floor(.4*RR(i-1));
    Pafp(i,3)=Rpeak(i)-floor(.1*RR(i-1));
    [~,I]=max(ecg(Pafp(i,1):Pafp(i,3))); Pafp(i,2)=I+Pafp(i,1)-1;
    
    Qafp(i,1)=Rpeak(i)-floor(.1*RR(i-1));    
    [~,I]=min(ecg(Qafp(i,1):Rpeak(i))); Qafp(i,2)=I+Qafp(i,1)-1;
    [~,I]=min(abs(ecg(Qafp(i,2):Rpeak(i)))); Qafp(i,3)=I+Qafp(i,2)-1;
%     Qafp(i,3)=Qafp(i,2)+floor(.01*RR(i-1));
%     Qafp(i,3)=Qafp(i,2)+round(.4*(Rpeak(i)-Qafp(i,2)));
    
    Safp(i,3)=Rpeak(i)+floor(.1*(RR(i)));    
    [~,I]=min(ecg(Rpeak(i):Safp(i,3))); Safp(i,2)=I+Rpeak(i)-1;
    [~,I]=min(abs(ecg(Rpeak(i):Safp(i,2)))); Safp(i,1)=I+Rpeak(i)-1;
%     Safp(i,1)=Safp(i,2)-floor(.01*RR(i));
    
    Rafp(i,1)=Qafp(i,3);  Rafp(i,3)=Safp(i,1);


    Tafp(i,3)=Rpeak(i)+floor(.6*RR(i));
    Tafp(i,1)=Rpeak(i)+floor(.1*RR(i));
    [~,I]=max(ecg(Tafp(i,1):Tafp(i,3))); Tafp(i,2)=I+Tafp(i,1)-1;
%     Tafp(i,2)=Rpeak(i)+floor(.4*RR(i));

end


Pafp(end,1)=Rpeak(end)-floor(.4*RR(end-1));
Pafp(end,3)=Rpeak(end)-floor(.1*RR(end-1));
[~,I]=max(ecg(Pafp(end,1):Pafp(end,3))); Pafp(end,2)=I+Pafp(end,1)-1;

Qafp(end,1)=Rpeak(end)-floor(.1*RR(end-1));    
[~,I]=min(ecg(Qafp(end,1):Rpeak(end))); Qafp(end,2)=I+Qafp(end,1)-1;
% % [~,I]=min(abs(ecg(Qafp(end,2):Rpeak(end)))); Qafp(end,3)=I+Qafp(end,2)-1;
% % Qafp(end,3)=Qafp(end,2)+floor(.01*RR(end-1));
Qafp(end,3)=Qafp(end,2)+round(.4*(Rpeak(end)-Qafp(end,2)));



Safp(1,3)=Rpeak(1)+floor(.1*(RR(1)));    
[~,I]=min(ecg(Rpeak(1):Safp(1,3))); Safp(1,2)=I+Rpeak(1)-1;
%     [~,I]=min(abs(ecg(Rpeak(1):Safp(1,2)))); Safp(1,1)=I+Rpeak(1)-1;
Safp(1,1)=Safp(1,2)-floor(.01*RR(1));

Rafp(end,1)=Qafp(end,3);  Rafp(1,3)=Safp(1,1);

Tafp(1,3)=Rpeak(1)+floor(.6*RR(1));
Tafp(1,1)=Rpeak(1)+floor(.1*RR(1));
[~,I]=max(ecg(Tafp(1,1):Tafp(1,3))); Tafp(1,2)=I+Tafp(1,1)-1;
% Tafp(1,2)=Rpeak(1)+floor(.4*RR(1));


end

