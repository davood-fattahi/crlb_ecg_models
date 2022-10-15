function [CRLB, FI]=GaussCrlbNumeric(varargin)

% Syntax: [CRLB, FI]=GaussCrlbNumeric(n, fs, p, varn, N );
%         [CRLB, FI]=GaussCrlbNumeric(t, p, varn, N);


if nargin==5
    n=varargin{1};
    fs=varargin{2};
    p=varargin{3};
    varn=varargin{4};
    N=varargin{5};
    t=n./fs;
elseif nargin==4
    t=varargin{1};
    p=varargin{2};
    varn=varargin{3};
    N=varargin{4};
end

p=p(:); t=t(:);
NumGaus=numel(p)./3;
if round(NumGaus)~=NumGaus %#ok<BDSCI>
    error('wrong dimension of p!')
end

dgdp=GausGrad(t,p); % gradient of gaussian respect to its parameters
    
FI=(N./varn)*dgdp*(dgdp');
CRLB=(varn./N).*inv(dgdp*(dgdp'));

end

