function afp =ecgfidapprox(ecg,fs,hr,varargin)

% This function appriximately finds peaks and borders of the input ecg beat
% waves.   

% afp =Ecgfidapprox(ecg,fs,hr,sm)

% Inputs:
% ecg: one ecg beat,
% fs: sampleing frequency
% hr: approximated heart rate in Hz
% sm: segmenting method, specified as 'separate', 'overlapped' and 'widest'.
% 'separate': the segments will have separate borders.
% 'overlapped': the segments are a little overlapped.
% 'widest': the segments are widely overlapped.


% Outputs:
% afp: the vector of approximated fiducal points


% Davood Fattahi, 1/12/2020
% fattahi.d@gmail.com




ecg=ecg(:);
[aR,R]=max(ecg);
[aQ,Q]=min(ecg(R-floor(.1*(fs/hr)):R)); Q=Q+R-floor(.1*(fs/hr))-1;
[aS,S]=min(ecg(R:R+floor(.1*(fs/hr)))); S=S+R-1;
[aP,P]=max(ecg(1:R-floor(.1*(fs/hr))));
[aT,T]=max(ecg(R+floor(.1*(fs/hr)):end-floor(.1*(fs/hr)))); T=T+R+floor(.1*(fs/hr))-1;
[aU,U]=max(ecg(T+floor(.1*(fs/hr)):end)); U=U+T+floor(.1*(fs/hr))-1;

if isequal(varargin{1},'separate')
    Ps=P-floor(.1*(fs/hr));
    Pe=P+floor(.1*(fs/hr));
    Qs=Pe+1;
    Qe=Q+ceil(abs(aQ./aR)*(R-Q));
    Rs=Qe+1;
    Ss=S-ceil(abs(aS./aR)*(S-R));
    Re=Ss-1;
    Se=S+ceil(abs(aS./aR)*(S-R));
    Ts=T-floor(.17*(fs/hr));
    Te=T+floor(.17*(fs/hr));
    Us=Te+1;
    Ue=size(ecg,1);
elseif isequal(varargin{1},'overlapped')  
    Ps=P-floor(.1*(fs/hr));
    Pe=P+floor(.1*(fs/hr));
    Qs=Q-floor(.05*(fs/hr));
    Qe=Q+ceil(abs(aQ./aR)*(R-Q));
    Rs=Qe;
    Ss=S-ceil(abs(aS./aR)*(S-R))-1;
    Re=Ss;
    Se=S+floor(.05*(fs/hr));
    Ts=T-floor(.17*(fs/hr));
    Te=T+floor(.17*(fs/hr));
    Us=Te;
    Ue=size(ecg,1);
elseif isequal(varargin{1},'widest')
    Ps=1;
    Pe=P+floor(.5*(Q-P));
    Qs=Pe;
    Qe=Q+floor(.5*(R-Q));
    Rs=Qe;
    Ss=S-floor(.5*(S-R));
    Re=Ss;
    Se=S+floor(.5*(T-S));
    Ts=Se;
    Te=T+floor(.17*(fs/hr));
    Us=Te;
    Ue=size(ecg,1);    
end
afp=[Ps Qs Rs Ss Ts Us ; P Q R S T U; Pe Qe Re Se Te Ue];


end

