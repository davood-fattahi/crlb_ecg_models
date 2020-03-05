function [ecgM,rrM]=ecgmean(ecg,fs,varargin)

% Syntax:
% [ecgM,rrM]=ecgmean(ecg,fs)
% [ecgM,rrM]=ecgmean(ecg,fs,'pan_tomp')
% [ecgM,rrM]=ecgmean(ecg,fs,'sameni',ff)

% Description:
% find the average ecg beat and average RR interval value, from the 1-d input
% ecg signal (common morphology of lead II). First, R peaks are detected, 
% then the beats are segmented by aligning the R peaks, while each beat 
% contains 35 percent of the previous RR interval and 65 percent of the 
% next one. No time warping is applied, and the segmented beats are 
% equalized in size using zeropadding. 
%
%
% Inputs:
% ecg: 1-d vector of ecg record
% fs: sampling frequency
% RPDM: R-peak detection method, specified as 'pan_tomp' for
% Pan-Tompin algorithm, or 'sameni' for Sameni algorithm.
% ff: esimated HR in Hz, only in RPDM = 'sameni'.  
% 
%
%
% Output:
% ecgM: the average ecg beat
% rrM: the average value of RR intervals
%
% Davood Fattahi, 08/01/2020
% fattahi.d@gmail.com




ecg=ecg-mean(ecg);
if(nargin==2 || strcmp(varargin{1},'pan_tomp'))
    [~,peakspos,~]=pan_tompkin(ecg,fs,0);
elseif(nargin==4 && strcmp(varargin{1},'sameni'))
    peaks = PeakDetection6(ecg,varargin{2},.3);
    peaks(end)=0;
    peakspos=floor(mean([find((peaks-[0 peaks(1:end-1)])==1);find((peaks-[0 peaks(1:end-1)])==-1)]));
else
      error('wrong inputs');
end



RRint=peakspos(2:end)-peakspos(1:end-1);
mx=max(RRint);
mn=min(RRint);
rrM=mean(RRint(:));
RRint=[RRint(1) RRint RRint(end)];
ecgB=zeros(size(peakspos,2),mx);
for i=1:size(peakspos,2)
    ss=peakspos(i)-floor(.35*RRint(i));
    if ss<1
        ss=1;
    end
    ee=peakspos(i)+floor(.65*RRint(i+1));
    if ee>size(ecg(:),1)
        ee=size(ecg(:),1);
    end
    s=ceil(0.35*mx) - peakspos(i)+ss;
    e=ceil(0.35*mx) -peakspos(i)+ee;
    ecgB(i,s:e)=ecg(ss:ee);
end

ecgM=mean(ecgB);

