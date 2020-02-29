function fp =EcgSegmentorApprox(ecg,fs,hr,varargin)

% This function appriximately segments the input ECG beat according to the R peak
% position.

% fp =EcgSegmentorApprox(ecg,fs,hr,ns)

% Inputs:
% ecg: one ecg beat,
% fs: sampleing frequency
% hr: approximate heart rate
% ns: number of segments, can be 3, 5, or 7.

% Outputs:
% fp: the vector of fiducal points


% developer: Davood Fattahi
% 1/12/2020




[~,R]=max(ecg);

if nargin==3 || varargin{1}==3 
    Ps=floor(R-.27*fs./hr);
    Pe=floor(R-0.0525*fs./hr);
    Qs=floor(R-0.0725*fs./hr);
    % Qe=floor(R-0.0125*fs./hr);
    % Rs=floor(R-0.0175*fs./hr);
    % Re=floor(R+0.0150*fs./hr);
    % Ss=floor(R+0.0125*fs./hr);
    Se=floor(R+.0625*fs./hr);
    STs=floor(R+.0525*fs./hr);
    % STe=floor(R+.17*fs./hr);
    % Ts=floor(R+.15*fs./hr);
    % Te=floor(R+.475*fs./hr);
    % Us=floor(R+.4*fs./hr);
    Ue=floor(R+.60*fs./hr);

    % fp=[Ps Qs Rs Ss STs Ts Us ; Pe Qe Re Se STe Te Ue];
    fp=[Ps Qs STs ; Pe Se Ue];
elseif nargin==4 && varargin{1}==5
    Ps=floor(R-.27*fs./hr);
    Pe=floor(R-0.0525*fs./hr);
    Qs=floor(R-0.0725*fs./hr);
    Qe=floor(R-0.0125*fs./hr);
    Rs=floor(R-0.0175*fs./hr);
    Re=floor(R+0.0150*fs./hr);
    Ss=floor(R+0.0125*fs./hr);
    Se=floor(R+.0625*fs./hr);
    STs=floor(R+.0525*fs./hr);
%     STe=floor(R+.17*fs./hr);
%     Ts=floor(R+.15*fs./hr);
%     Te=floor(R+.475*fs./hr);
%     Us=floor(R+.4*fs./hr);
    Ue=floor(R+.60*fs./hr);

    fp=[Ps Qs Rs Ss STs ; Pe Qe Re Se Ue];
    
elseif nargin==4 && varargin{1}==7
    Ps=floor(R-.27*fs./hr);
    Pe=floor(R-0.0525*fs./hr);
    Qs=floor(R-0.0725*fs./hr);
    Qe=floor(R-0.0125*fs./hr);
    Rs=floor(R-0.0175*fs./hr);
    Re=floor(R+0.0150*fs./hr);
    Ss=floor(R+0.0125*fs./hr);
    Se=floor(R+.0625*fs./hr);
    STs=floor(R+.0525*fs./hr);
    STe=floor(R+.17*fs./hr);
    Ts=floor(R+.15*fs./hr);
    Te=floor(R+.475*fs./hr);
    Us=floor(R+.4*fs./hr);
    Ue=floor(R+.60*fs./hr);

    fp=[Ps Qs Rs Ss STs Ts Us ; Pe Qe Re Se STe Te Ue];
    
else 
    error('wrong input arguments')
end

