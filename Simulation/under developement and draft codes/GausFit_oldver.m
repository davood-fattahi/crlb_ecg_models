function [p, varargout]=GausFit(t,x,p0, lb, ub, varargin)
% Fit sum of Gaussian functions on the data.
% 
% Syntax:
% [p, varargout]=GausFit(t,x,p0, lb, ub)
% [p, varargout]=GausFit(t,x,p0, lb, ub, options)
% [p, varargout]=GausFit(t,x,p0, lb, ub, PrMu, PrCov, PstCov)
% [p, varargout]=GausFit(t,x,p0, lb, ub, PrMu, PrCov, PstCov, options)
% 
% Inputs:
% t: time samples
% x: input data
% p0: parameters' initial values. The first, second, and third rows contain
%   respectively Gaussians' amplitude (a), Gaussians' width (b), and
%   Gaussians' center (c). So the i-th column contains the i-th Gaussian
%   parameters, by the order of [ai ; bi ; ci]. It also can be in the
%   vector form, like [a1 b1 c1 a2 b2 c2 ...]'.
% lb: lower bound for parameters
% ub: upper bound for parameters
% options: options structure for optimization problem. (Note: It should be in a
%   simple strutre format, NOT the 'optimoptions' format commonly used in matlab
%   optimization toolbox. However, the fields are same as optimoptions.)
% PrMu: mean of prior distribution
% PrCov: covarinace matrix of prior information
% PstCov: covariance or variane of noise
% 
% Output:
% p: the estimated parameters of Gaussians, with the same format as p0;
%
% Davood Fattahi, 3/31/2021, last edition 4/27/2021.
% fattahi.d@gmail.com

p0=p0(:); t=t(:); x=x(:); lb=lb(:); ub=ub(:);
N=length(x); NumGaus=length(p0)/3;
if round(NumGaus)~=NumGaus
    error('Wrong number of parameters in p0!')
end

if length(varargin)<=1 % ML
    ML=true; Bys=false;
elseif length(varargin)>=3 % Bys
    ML=false; Bys=true;
    PrMu=varargin{1}; PrMu=PrMu(:);
    PrCov=varargin{2};
    PstCov=varargin{3};
    if isscalar(PstCov)
        PstCov=PstCov.*eye(N);
    end
end
if nargin==5 || nargin==8
    options0=struct;
elseif nargin==6 || nargin==9
    options0=varargin{end};
end


%%% Gaussfit

if ML % ML
    fun=@(p)(GausVal(t,p)-x);
elseif Bys % Bys
    fun=@(p)([sqrtm(PstCov)\(GausVal(t,p)-x) ; sqrtm(PrCov)\(p-PrMu)]);
end
options = optimoptions(@lsqnonlin,'Display','off');
options=optionsmerge(options, options0);
p = lsqnonlin(fun,p0,lb,ub,options);


    


function options=optionsmerge(options1, options2)
%%% \\\\\\\\\\\\\\\\\\\\\\
%%% merging options2 into options1.
FN=fieldnames(options2);
for i=1:length(FN)
    eval(['options1.' FN{i} '=options2.' FN{i} ';'])
end
options=options1;

