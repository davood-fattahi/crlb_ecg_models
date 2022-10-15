function  [signal, newPeaks]=ECGResembler(tm, pmean, pstd, rrmean, rrstd, varargin)

% signal=ECGResembler(tm, pmean, pstd, rrmean, rrstd, phase0)
% signal=ECGResembler(tm, pmean, pstd, rrmean, rrstd)

% The pmean, rrmean, and rrstd fields can be empty, then they are calculated from the input sample ecg:
% signal=ECGResembler(tm, pmean, pstd, rrmean, rrstd, phase0, ecg, ecgFs)
% signal=ECGResembler(tm, pmean, pstd, rrmean, rrstd, phase0, ecg, ecgFs,ecgHr)
% signal=ECGResembler(tm, pmean, pstd, rrmean, rrstd, phase0, ecg, ecgFs, ecgHr, phasebins)
% Note: pstd cannot be empty.

% pstd: parameters relative standard deviation around their mean.
% rrstd: rr interval relative standard deviation around its mean.

%% calculating the resembling parameters

if nargin==5
    phase0=0;
elseif nargin==6
    phase0=varargin{1};
elseif nargin==8
    phase0=varargin{1};
    ecg=varargin{2};
    ecgFs=varargin{3};
    ecgHr=72/60/ecgFs;
    phasebins=300;
elseif nargin == 9
    phase0=varargin{1};
    ecg=varargin{2};
    ecgFs=varargin{3};
    ecgHr=varargin{4};
    phasebins=300;
elseif nargin == 10
    phase0=varargin{1};
    ecg=varargin{2};
    ecgFs=varargin{3};
    ecgHr=varargin{4};
    phasebins=varargin{5};    
else
    error('the number of input arguments are not consistent!');
end

tm=tm(:);
if nargin>=8
%     % baseline wandering removal,
%     ecg=ecg-(BaseLine1(BaseLine1(ecg', round(w1*ecgFs), 'md'), round(w2*ecgFs), 'mn'))';

    % extract the R peaks using peak detector
    Peaks = PeakDetection(ecg,ecgHr); Rp=find(Peaks);
    
    if isempty(rrmean) 
        rrmean=mean(Rp(2:end)-Rp(1:end-1),'all')/ecgFs;
    end
    if isempty(rrstd)
        rrstd=std(Rp(2:end)-Rp(1:end-1),0,'all')/ecgFs/rrmean;
    end
    if isempty(pmean)
        % phase calculation
        [phase, ~] = PhaseCalculation(Peaks);     % phase calculation

        % mean ecg
        [ECGmean,~,meanphase] = MeanECGExtraction(ecg,phase,phasebins,1); % mean ECG extraction

        % fit 7 gaussians
        p0=[.03 -0.03 -.1  1  -.1  .15  .15; .2 .2 .05 .05 .05 .3 .3; -1.5 -1.3 -.2  0   .2  1.8 2.2]; 
        lb=[0   -inf  -inf 0  -inf -inf  -inf; 0  0  0   0   0   0  0; -2   -2   -.5  -.1  0  .5  1]; 
        ub=[inf  inf  0    inf  0    inf  inf; 1  1 .15 .1  .1  1.2 1.2;  -.5 -.5  0    .1  .5  2.5 2.5];
        pmean=GausFit(meanphase,ECGmean, p0, lb, ub, struct('SpecifyObjectiveGradient',true));
% 
%         figure; plot(meanphase,ECGmean); hold on
%         plot(meanphase,GausVal(meanphase,pmean))
    end
end


%% resembling
pmean=reshape(pmean,3,[]);

if isscalar(pstd)
    pstd=pstd.*pmean;
end

newPeaks=zeros(size(tm));
signal=zeros(size(tm));
fs=round((length(tm)-1)./(tm(end)-tm(1)));

RRmean=floor(rrmean*fs); 
RRstd=floor(rrstd*rrmean*fs); 

NumBeats=floor(2*(range(tm)./rrmean));
start = floor((mod(phase0 + pi, 2*pi)/(2*pi))*RRmean);
Rp=cumsum([start; floor(RRmean+RRstd*randn(NumBeats,1))]);
Rp(Rp>length(tm))=[];
newPeaks(Rp)=true;

[newphase, ~] = PhaseCalculation(newPeaks);     % phase calculation

p=repmat(pmean,1,1,NumBeats)+pstd.*randn(size(pmean,1),size(pmean,2), NumBeats); % parameters deviation
for i=1:length(Rp)
    p(2,(p(2,:,i)<=0),i)=pmean(2,(p(2,:,i)<=0)); % check if b values are negative
    if i==1
        [signal(1:Rp(i)-1), ~]= SingleChannelECGGenerator(newphase(1:Rp(i)-1),0,p(1,:,i),p(2,:,i),p(3,:,i));        
    else
        [signal(Rp(i-1):Rp(i)-1), ~]= SingleChannelECGGenerator(newphase(Rp(i-1):Rp(i)-1),0,p(1,:,i),p(2,:,i),p(3,:,i));
    end
end


