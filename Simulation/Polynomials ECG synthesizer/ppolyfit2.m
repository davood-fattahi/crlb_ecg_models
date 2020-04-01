function [P, pcs, pct] =ppolyfit2(t,x,PolyOrder,varargin)

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
% [P, pcs, pct]=ppolyfit(t,x,PolyOrder,pcs,ext);
% fit polynomials on each piece of the function x having values at time 
% vector of t. The pieces  start and end sample indexes are given by pcs
% vector, but in the optimization each piece is extended by length of ext to
% reduce roughness of the borders.
% 
% [P, pcs, pct]=ppolyfit(t,x,PolyOrder,ws,nvrlp,ext);
% fit polynomials on each piece of the function x having values at time 
% vector t. Each piece's size is given by ws in sample, with overlap 
% samples determined in nvrlp.
%
% Input: 
% t: the vestor of time
% x: the values vector
% PolyOrder: order of polynomials, can be a scalar or a vector with the 
% size equal to number of the pieces.
% pcs: an N by 2 matrix in which N is the number of pieces and each row
% contains the start and end sample indexes of each piece.
% nvrlp: overlap size of pieces. 
% ext: extention size in sample, a scalar.
% 
% Outputs:
% P: the polynoials coefficients in cell format.
% pcs: the pieces start and end samples' indexes.
% pct: the pieces start and end time.
% 


xsz=size(x(:),1); % signal size in samples
if nargin==3
    ext=0;
    ws=floor(xsz./5);
    nvrlp=floor(ws./2);
    np=ceil((xsz-ws)./(ws-nvrlp))+1;
    for i=1:np
    pcs(i,:)=[(i-1)*(ws-nvrlp)+1 (i-1)*(ws-nvrlp)+ws];
    end
elseif nargin==4
    ext=0;
    pcs=varargin{1};
    np=size(pcs,1);
elseif nargin==5
    ext=varargin{2};
    pcs=varargin{1};
    np=size(pcs,1);    
elseif nargin==6
    ext=varargin{3};
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
    ws=pcs(i,2)-pcs(i,1)+1;
    if ext<1 && ext>0
        ext=floor(ext.*ws);
    end
    pci=pcs(i,1)-ext:pcs(i,2)+ext;
    pci(pci<1|pci>xsz)=[];
    xx=x(pci); xx=xx(:);
    tt=t(pci); tt=tt(:); pct(i,:)=[tt(1) tt(end)];

    
    %%% range normalization
    beta=2/(tt(end)-tt(1));
    gamma=-(2*tt(1)/t(end))-1;
    ttt=tt.*beta+gamma;
    
    %%% amplitude normalization
    Beta=2./(max(xx)-min(xx));
    Gamma=-1;
    xxx=xx.*Beta+Gamma;
    
    %%% polyfit
    [p,~]=polyfit(ttt,xxx,PolyOrder(i));
    
    %%% range denormalization
    coefs=flip(p);
    C=zeros(PolyOrder(i)+1);
    for j=0:PolyOrder(i)
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

