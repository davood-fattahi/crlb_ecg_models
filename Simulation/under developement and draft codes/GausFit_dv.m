function [p, varargout]=GausFit2(t,x,p0, lb, ub, varargin)

% demo code for fittng Gaussuans + line ...




p0=p0(:); t=t(:); x=x(:); lb=lb(:); ub=ub(:);
N=length(x); NumGaus=length(p0)/5;

if length(varargin)<=1 % ML
    ML=true; Bys=false;
elseif length(varargin)>=3
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
g='@(p)';
for i=1:NumGaus
    g=([g '+(p('  num2str((i-1)*5+1) ')*(t-p('  num2str((i-1)*5+2) '))+p(' num2str((i-1)*5+3) ').*exp(-((t-p(' num2str((i-1)*5+5) ')).^2)./(2*p(' num2str((i-1)*5+4) ')^2)))']);
end
g=eval(g);

if ML % ML
    fun=@(p)(g(p)-x);
    options = optimoptions(@lsqnonlin,'Display','off');
    options=optionsmerge(options, options0);
    p = lsqnonlin(fun,p0,lb,ub,options);
elseif Bys % Bys
    fun=@(p)((g(p)-x)'*((PstCov)\(g(p)-x))+(p-PrMu)'*(PrCov\(p-PrMu)));
    if isempty(lb) && isempty(ub)
        options = optimoptions(@fminunc,'Display','off');
        options=optionsmerge(options, options0);
        p=fminunc(fun,p0,options);
    else
        options = optimoptions(@fmincon,'Display','off');
        options=optionsmerge(options, options0);
        p = fmincon(fun,p0,[],[],[],[],lb,ub,[],options);
    end
end


    


function options=optionsmerge(options1, options2)
FN=fieldnames(options2);
for i=1:length(FN)
    eval(['options1.' FN{i} '=options2.' FN{i} ';'])
end
options=options1;

