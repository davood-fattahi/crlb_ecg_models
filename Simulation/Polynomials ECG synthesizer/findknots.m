function k=findknots(ecg,depth,varargin)
% 
% Syntax:
% k=ecgknots(ecg,depth,tr)
% 
% Description:
% Finding the suitable knots along the ecg signal for spline or polynomial
% fitting. At the first level (depth=0), the first and last samples are
% set as initial knots, and for each next levels, the point with maximum 
% distance from the reference line (the line passing through the knots) are 
% considered as a new knot.
% 
% Inputs:
%   ecg: A beat of ecg signal.
%   depth: the number of levels to find the knots within the segments.
%   tr: A factor of variance as the minimum teshold of diatance for which 
%       the corresponding knot will be kept. For example; if tr=2 it means
%       the found new knot is discarded if its distance from the refrence
%       line is smaller than 2*variance of the segment of signal between 
%       the two neighbor knots. The default value is 1.
% 
% Output:
% k: Vector of the found knots
% 
% 
% Davood Fattahi, 10/02/2020;
% fattahi.d@gmail.com
% 
% 

if nargin==3
    tr=varargin{1};
elseif nargin==2
    tr=0;
else
    error('wrong nuber of inputs')
end

ecg=ecg(:)';
k=nan(2^depth,2);
kk=k;
k(1,1)=1;
k(1,2)=size(ecg(:),1);
for j=1:depth
for i=1:2^(j-1)
    tt=k(i,1):k(i,2);
    s=ecg(tt);
    line=s(1)+(1:size(s(:),1)).*((s(end)-s(1))/size(s(:),1));
    [M,I]=max(abs(s-line));
    if M>tr*(std(s))
    kk(2*i-1,1)=k(i,1);
    kk(2*i-1,2)=tt(I);
    kk(2*i,1)=tt(I);
    kk(2*i,2)=k(i,2);
    else
    kk(2*i-1,1)=k(i,1);
    kk(2*i-1,2)=k(i,2);
    kk(2*i,1)=k(i,1);
    kk(2*i,2)=k(i,2);        
    end
end
k=kk;
kk=nan(2^depth,2);
end
k=k(:);
k(isnan(k))=[];
k=unique(k);

