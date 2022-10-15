function [coefs, varargout]=PolyFit(t,x,PolyOrder,scflag,varargin)
%
% Syntax:
% [coefs, sccoefs]=PolyFit(t,x,PolyOrder,scflag, PrMu, PrCov, PstCov)
% coefs=PolyFit(t,x,PolyOrder,scflag, PrMu, PrCov, PstCov)
% Bayesian estimation of polynomials coeficients, using Tikhonove
% regularized linear least squares.
% 
% [coefs, sccoefs]=PolyFit(t,x,PolyOrder,scflag)
% coefs=PolyFit(t,x,PolyOrder,scflag)
% ML estimation of polynomials coeficients, using linear least squres.
% 
% Inputs:
% t: vector of time samples,
% x: values for each time sample
% PolyOrder: order of polynomilas
% scflag: if true, t and x ranges are scaled and centeralized to [-1 1] 
%         befor estimating the coefs, the final delivered coefs will be
%         de-scaled and de-centeralized. 
% PrMu: mean of coefs,
% PrCov: covariance matrix of coefs,
% PstCov: variance of noise 
% 
% Outputs:
% coefs: polynomials coefficients.
% sccoefs: polynomials coefficients for scaled and centeralized x and t.

% davood fattahi, fattahi.d@gmail.com
% 3 Feb 2021

x=x(:); t=t(:);
N=length(x);

if scflag
    %%% range normalization
    beta=2/(t(end)-t(1));
    gamma=-(t(1)+t(end))/(t(end)-t(1));
    t=t.*beta+gamma;

    %%% amplitude normalization
    Beta=2./(max(x)-min(x));
    Gamma=-(max(x)+min(x))/(max(x)-min(x));
    x=x.*Beta+Gamma;
end

if isempty(varargin) % ML
    ML=true; Bys=false;
else
    ML=false; Bys=true;
    PrMu=varargin{1}; PrMu=PrMu(:);
    PrCov=varargin{2};
    PstCov=varargin{3};
    if isscalar(PstCov)
        PstCov=PstCov.*eye(N);
    end
end


%%% polyfit
H=zeros(PolyOrder+1,N);
for i=0:PolyOrder
    H(i+1,:)=t.^i;
end
if ML % ML
    p=(H*H')\(H*x);
elseif Bys % Bys
    A=H*(PstCov\H')+inv(PrCov);
    B=H*(PstCov\(x-H'*PrMu));
    p=PrMu+A\B;
end

if scflag
    %%% range denormalization
    C=zeros(PolyOrder+1);
    for j=0:PolyOrder
        for k=0:j
            C(j+1,k+1)=nchoosek(j,k)*(beta^k)*(gamma^(j-k))*p(j+1);
        end
    end
    coefs=sum(C,1);

    %%% amplitude denormalization
    coefs(1)=coefs(1)-Gamma;
    coefs=coefs./Beta;
    sccoefs=p;
else
    coefs=p;
    sccoefs=[];
end


if nargout>1
    varargout{1}=sccoefs;
end

    