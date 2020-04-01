function [S, SS]=ppolyval(P,pct,varargin)

% calculate piecewise polynomials values.
% 
% [S, SS]=ppolyval(P,pct,sr)
% 
% Inputs:
% P: a cell format containing polynomials coefficients.
% pct: each of the pieces time vector. It can be an N by 2 vector in which 
% the first column is the start times and the second column is the end
% time. In this case, the sample rate (sr) in necessary in inputs. In 
% another case pct is a matrix in which each row demonstrate the time 
% vector of each piece. 
% sr: sampling rate.
%
% Outputs:
% S: the vector of all values averages over overlapped intervals.
% SS: the matrix of values in which each row is correspounding to each
% piece.
% 
% Davood Fattahi, 01/03/2020
% fatahi.d@gmail.com
%
%

if nargin==3
    sr=varargin{1};
    ssz=ceil((pct(end)-pct(1)).*sr);    
elseif nargin==2
    ssz=size(unique(pct(:)),1);
end
SS=nan(size(pct,1),ssz);

for i=1:size(pct,1)
    if size(pct,2)==2
        tt=pct(i,1):1/sr:pct(i,end); tt=tt(:);
    elseif size(pct,2)>2
        tt=pct(i,:); tt=tt(:);
    else
        error 'wrong dimension of pcs!'
    end
    pci=(floor(pct(i,1)*sr)+1:floor(pct(i,1)*sr)+size(tt,1))-floor(pct(1,1)*sr);
    S=polyval(P{i,1},tt);
    SS(i,pci)=S;
end
S=mean(SS,1,'omitnan');
