function [ecgM,rrM]=ecgmean(ecg,fs,varargin)

ecg=ecg-mean(ecg);
if(nargin==2 || strcmp(varargin{1},'pan_tomp'))
    [~,peakspos,~]=pan_tompkin(ecg,fs,0);
elseif(nargin==4 && strcmp(varargin{1},'sameni'))
    peaks = PeakDetection6(ecg,varargin{2},max(ecg)./5);
    peaks(end)=0;
    peakspos=floor(mean([find((peaks-[0 peaks(1:end-1)])==1);find((peaks-[0 peaks(1:end-1)])==-1)]));
else
      error('wrong inputs');
end



RRint=peakspos(2:end)-peakspos(1:end-1);
mx=max(RRint);
mn=min(RRint);
rrM=mean(RRint(:));
ecgB=zeros(size(peakspos,2),mx);
for i=2:size(peakspos,2)-1
    ss=peakspos(i)-floor(.33*RRint(i-1));
    ee=peakspos(i)+floor(.67*RRint(i));
    s=ceil(0.33*mx) - floor(.33*RRint(i-1));
    e=ceil(0.33*mx) + floor(.67*RRint(i));
    ecgB(i,s:e)=ecg(ss:ee);
end
ecgM=mean(ecgB);  

