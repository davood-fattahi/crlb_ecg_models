function Y = ECGsplinefit(ECG,type,order,varargin)

%%% range noramalization
t=normalize((1:size(ECG(:),1)),'range').*2-1;


if isstring(varargin{1})
    knots=findknots(ECG,varargin{1});
elseif isnumeric(varargin{1})
    knots=varargin{1};
end

Y=interp1(t,ECG,t)

%%% range denormalization
