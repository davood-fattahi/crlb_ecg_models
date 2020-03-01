function [PCoef, Y] = ECGsplinefit(ECG,type,order,varargin)
% 
% Syntax:
% [PC, Y] = ECGsplinefit(ECG,type,order)
% [PC, Y] = ECGsplinefit(ECG,type,order,xrange)
% [PC, Y] = ECGsplinefit(ECG,type,order,xrange,knotmethod)
% [PC, Y] = ECGsplinefit(ECG,type,order,xrange,knotmethod,...)
% 
% Description:
% 
% 
% Inputs:
% 
% 
% Outputs:
% 
% 
% 
% 
% Davood Fattahi, 20/02/2020
% fattahi.d@gmail.com


%%

%%% range noramalization
t=normalize((1:size(ECG(:),1)),'range',[-1,1]);


if isstring(varargin{1})
    knots=findknots(ECG,varargin{1});
elseif isnumeric(varargin{1})
    knots=varargin{1};
end

Y=interp1(t,ECG,t)

%%% range denormalization

    
for i=0:order
    for j=0:i
        C(i+1,j+1)=nchoosek(i,j)*(beta^j)*(gamma^(i-j))*coefs(i+1);
    end
end
coefs=sum(C,1);





