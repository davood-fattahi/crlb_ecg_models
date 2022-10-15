function p=GausFitCaruana(t,x, varargin)
% Fit a Gaussian functions on the data, using Caruana's method.
% 
% Syntax:
% [p, varargout]=GausFitPolyAppr(t,x)
% [p, varargout]=GausFitPolyAppr(t,x, PrMu, PrCov, PstCov)
% 
% Inputs:
% t: time samples
% x: input data
% PrMu: mean of prior distribution
% PrCov: covarinace matrix of prior information
% PstCov: covariance or variane of noise
% 
% Output:
% p: the estimated parameters of Gaussians.
%
% Davood Fattahi, 4/20/2021
% fattahi.d@gmail.com

[~, I]=max(abs(x));

if x(I)<0
    nflag=1;
    x=-x;
else
    nflag=0;
end

x(x<=eps)=eps;
t=t(:); x=log(3+x(:));

if nargin==2 % ML
    coefs=PolyFit(t,x,2,1);
elseif nargin==5 % Bys
    PrMu=varargin{1}; PrMu=PrMu(:);
    PrCov=varargin{2};
    PstCov=varargin{3};
    if isscalar(PstCov)
        PstCov=PstCov.*eye(N);
    end
    coefs=PolyFit(t,x,2,1,PrMu,PrCov,PstCov);
end
p=nan(3,1);
p(2)=(sqrt(-1/(2*coefs(3))));
p(3)=-coefs(2)/(2*coefs(3));
p(1)=exp(coefs(1)-(coefs(2)^2)/(4*coefs(3)));

if nflag
    p(1)=-p(1);
end





