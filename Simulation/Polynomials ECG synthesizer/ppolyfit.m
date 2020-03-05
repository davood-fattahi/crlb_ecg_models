function [P, pcs] =ppolyfit(t,x,PolyOrder,varargin)

% piecewise polynomials fitter
% 
% Syntax:
% P=ppolyfit(t,x,PolyOrder);
% P=ppolyfit(t,x,PolyOrder,pcs);
% P=ppolyfit(t,x,PolyOrder,ws,nvrlp);


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
    tt=t(pci); tt=tt(:);
    
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

