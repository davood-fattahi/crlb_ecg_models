function [P, pcs, pct] =ppolyfit(t,x,PolyOrder,varargin)

% Description:
% piecewise polynomials fitter
% 
% Syntax:
% [P, pcs, pct]=ppolyfit(t,x,PolyOrder);
% fit polynomials on each piece of the function x having values at time 
% vector t. Each piece's size is 20 percent of x length in sample, with 50
% percent overlap. PolyOrder determines the order of polynomials.
%  
% [P, pcs, pct]=ppolyfit(t,x,PolyOrder,pcs);
% fit polynomials on each piece of the function x having values at time 
% vector of t. The pieces  start and end sample indexes are given by pcs
% vector.
% 
% [P, pcs, pct]=ppolyfit(t,x,PolyOrder,ws,nvrlp);
% fit polynomials on each piece of the function x having values at time 
% vector t. Each piece's size is given by ws in sample, with overlap 
% samples determined in nvrlp.
%
% Input: 
% t: the vestor of time
% x: the values vector
% PolyOrder: order of polynomials, can be a scalar or a vector with the 
% size equal to number of the pieces.
% pcs: a N by 2 vector in which N is the number of pieces and each row
% contains the start and end sample indexes of each piece.
% nvrlp: overlap size of pieces. 
% ext: extention size in sample, a scalar.
% 
% Outputs:
% P: the polynoials coefficients in cell format.
% pcs: the pieces start and end samples' indexes.
% pct: the pieces start and end time.
% 
%
% Davood Fattahi, 01/03/2020
% fattahi.d@gmail.com


xsz=size(x(:),1); % signal size in samples
if nargin==3
    ws=floor(xsz./5);
    nvrlp=floor(ws./2);
    np=ceil((xsz-ws)./(ws-nvrlp))+1;
    for i=1:np
    pcs(i,:)=[(i-1)*(ws-nvrlp)+1 (i-1)*(ws-nvrlp)+ws];
    end
elseif nargin==4
    pcs=varargin{1};
    np=size(pcs,1);
elseif nargin==5
    ws=varargin{1};
    nvrlp=varargin{2};
    np=ceil((xsz-ws)./(ws-nvrlp))+1;
    for i=1:np
    pcs(i,:)=[(i-1)*(ws-nvrlp)+1 (i-1)*(ws-nvrlp)+ws];
    end
else
     error 'wrong inputs!'
end


PolyOrder=PolyOrder(:);
if size(PolyOrder,1)==1
    PolyOrder=PolyOrder.*ones(np,1);
elseif size(PolyOrder,1)~=np
    error 'wrong size of PolyOrder vector'
end


%%% correcting pieces
pcs(pcs<1)=1;
pcs(pcs>xsz)=xsz;

P=cell(np,1);
for i=1:np
    pci=pcs(i,1):pcs(i,2);

    xx=x(pci); xx=xx(:);
    tt=t(pci); tt=tt(:); pct(i,:)=[tt(1) tt(end)];
    
    %%% range normalization
    beta=2/(tt(end)-tt(1));
    gamma=-1;
    ttt=tt.*beta+gamma;
    
    %%% amplitude normalization
    Beta=2./(max(xx)-min(xx));
    Gamma=-1;
    xxx=xx.*Beta+Gamma;
    
    %%% polyfit
    [p,~]=polyfit(ttt,xxx,PolyOrder(i));
    
    %%% range denormalization
    coefs=flip(p);
    for j=0:PolyOrder
        for k=0:j
            C(j+1,k+1)=nchoosek(j,k)*(beta^k)*(gamma^(j-k))*coefs(j+1);
        end
    end
    coefs=sum(C,1);
    
    %%% amplitude denormalization
    coefs(1)=coefs(1)-Gamma;
    coefs=coefs./Beta;
    
    P{i,1}=flip(coefs);

end

