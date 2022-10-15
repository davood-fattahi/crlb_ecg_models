function [p]=GausFitAuto(t,x,NumGaus, varargin)
% Fit sum of Gaussian functions on the data.
% 
% Syntax:
% [p, varargout]=GausFitAuto(t,x,NumGaus)
% [p, varargout]=GausFitAuto(t,x,NumGaus, options)
% [p, varargout]=GausFitAuto(t,x,NumGaus, PrMu, PrCov, PstVar)
% [p, varargout]=GausFitAuto(t,x,NumGaus, PrMu, PrCov, PstVar, options)
% 
% Inputs:
% t: time samples
% x: input data
% p0: parameters' initial values. The first, second, and third rows contain
%   respectively Gaussians' amplitude (a), Gaussians' width (b), and
%   Gaussians' center (c). So the i-th column contains the i-th Gaussian
%   parameters, by the order of [ai ; bi ; ci]. It also can be in the
%   vector form, like [a1 b1 c1 a2 b2 c2 ...]'.
% NumGaus: number of gaussians.
% options: options structure for optimization problem. (Note: It should be in a
%   simple strutre format, NOT the 'optimoptions' format commonly used in matlab
%   optimization toolbox. However, the fields are same as optimoptions.)
% PrMu: mean of prior distribution
% PrCov: covarinace matrix of prior information
% PstVar: variane of the noise
% 
% Output:
% p: the estimated parameters of Gaussians in the vector form, like [a1 b1 c1 a2 b2 c2 ...]'.
%
% Davood Fattahi, 5/5/2021.
% fattahi.d@gmail.com

t=t(:); x=x(:);


if length(varargin)<=1 % ML
    ML=true; Bys=false;
elseif length(varargin)>=3 % Bys
    ML=false; Bys=true;
    PrMu=varargin{1}; PrMu=PrMu(:);
    PrCov=varargin{2};
    PstVar=varargin{3};
end
if nargin==3 || nargin==6
    options0=struct;
elseif nargin==4 || nargin==7
    options0=varargin{end};
end

%% initial point search
xx=x; p0=zeros(3,NumGaus);
for i=1:NumGaus
    [~,I]=max(abs(xx)); p0(1,i)=xx(I); p0(3,i)=t(I);
    p0(2,i) = lsqnonlin(@(bb)(MlFun(t,[p0(1,i); bb; p0(3,i)],xx)),range(t)/20,0,range(t)/5);
    xx=xx-GausVal(t,p0(:,i));
    zi=I-4*floor(p0(2,i)/(t(2)-t(1))):I+4*floor(p0(2,i)/(t(2)-t(1))); zi(zi<1 | zi>length(xx))=[];
    xx(zi)=0;
end


%%% Gaussfit

if ML % ML
    fun=@(p)(MlFun(t,p,x));
elseif Bys % Bys
    fun=@(p)(BysFun(t,p,x,PstVar, PrCov, PrMu));
end
options = optimoptions(@lsqnonlin,'Display','off');
options=optionsmerge(options, options0);
p = lsqnonlin(fun,p0,[],[],options); 
p=p(:);


end

%% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%%%% functions

function [F, J]=MlFun(t,p,x)
    F=GausVal(t,p)-x;
    if nargout>1
        J=GausGrad(t,p)';
    end
end

function [F, J]=BysFun(t,p,x,PstVar, PrCov, PrMu)
    F=([(GausVal(t,p)-x)./sqrt(PstVar) ; sqrtm(PrCov)\(p-PrMu)]);
    if nargout>1
        J=[GausGrad(t,p)'./sqrt(PstVar) ; inv(sqrtm(PrCov))];
    end
end

function options=optionsmerge(options1, options2)
%%% \\\\\\\\\\\\\\\\\\\\\\
%%% merging options2 into options1.
FN=fieldnames(options2);
for i=1:length(FN)
    eval(['options1.' FN{i} '=options2.' FN{i} ';'])
end
options=options1;
end
