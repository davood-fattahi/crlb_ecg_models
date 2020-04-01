function varargout=polydeviate(pp,beta,gamma,amp,varargin)

% deviate the polynomial coefficients according to this parameter changing:
% amp*(beta*x+gamma)--->x
% 
% 
% [coefs,tt] = polydeviate(pp,beta,gamma,amp,t);
% coefs = polydeviate(pp,beta,gamma,amp);



PO=size(pp,2);
coefs=flip(pp,2);
for i=1:size(coefs,1)
C=zeros(PO);
for j=0:PO-1
    for k=0:j
        C(j+1,k+1)=nchoosek(j,k)*(beta^k)*(gamma^(j-k))*coefs(i,j+1);
    end
end
coefs(i,:)=sum(C,1);
end

coefs=flip(coefs,2);
coefs=coefs*amp;

varargout{1}=coefs;

if nargin==5
    varargout{2}=(varargin{1}-gamma)/beta;
end
end
