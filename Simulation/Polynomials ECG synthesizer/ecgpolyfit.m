function [P, AF] =ecgpolyfit(ecg,hr,fs,PolyOrder,varargin)

% 
% 
% [P, AF] =ecgpolyfit(ecg,hr,fs,PolyOrder)
%
% [P, AF] =ecgpolyfit(ecg,hr,fs,PolyOrder,Pcs)
%     Pcs = {'P,Q,R,S,T','P,Q,R,S,T,U','P,QRS,T','P,QRS,T,U'}
% 
% [P, AF] =ecgpolyfit(ecg,hr,fs,PolyOrder,Pcs,FAM)
% 
% [P, AF] =ecgpolyfit(ecg,hr,fs,PolyOrder,Pcs,FAM,PFM)
% 
% [P, AF] =ecgpolyfit(ecg,hr,fs,PolyOrder,'MovingWindow')
% 
% [P, AF] =ecgpolyfit(ecg,hr,fs,PolyOrder,'MovingWindow',ws)
% 
% [P, AF] =ecgpolyfit(ecg,hr,fs,PolyOrder,'MovingWindow',ws,ovrlp)
% 
% 
% 

ss=size(ecg(:),1); % ecg size in samples

%%% setting default heart rate
if isempty(hr)
    hr=fs./ss;
end
%%%

if nargin<4
    error 'Not enough inputs!'
elseif nargin==4
    Pcs='MovingWindow';
    ws=floor(.2*ss);
    nvrlp=floor(.5*ws);
elseif nargin==5
    Pcs=varargin{1};
    switch Pcs
        case {'P,Q,R,S,T','P,Q,R,S,T,U','P,QRS,T','P,QRS,T,U'}
            FAM='overlapped';
            PFM=[];
        case 'MovingWindow'
            ws=floor(.2*ss);
            nvrlp=floor(.5*ws);
    end
elseif nargin==6
    Pcs=varargin{1};
    switch Pcs
        case {'P,Q,R,S,T','P,Q,R,S,T,U','P,QRS,T','P,QRS,T,U'}
            FAM=varargin{2};
            PFM=[];
        case 'MovingWindow'
            error 'Not enough inputs! both the window size and overlap ratio is required.'
    end
elseif nargin==7
    Pcs=varargin{1};
    switch Pcs
        case {'P,Q,R,S,T','P,Q,R,S,T,U','P,QRS,T','P,QRS,T,U'}
            FAM=varargin{2};
            PFM=varargin{3};
        case 'MovingWindow'
            ws=varargin{2};
            nvrlp=floor(varargin{3}.*ws);
    end 
else
    error 'too many inputs!'
        
end




if isequal(Pcs,'MovingWindow')
    np=ceil((ss-ws)./(ws-nvrlp))+1;
    figure
    plot(ecg)
    hold on
    for i=1:np
        t=(i-1)*(ws-nvrlp)+1:(i-1)*(ws-nvrlp)+ws;
        if t(end)>ss
            t=(i-1)*(ws-nvrlp)+1:ss;
        end
        s=ecg(t);
        s=s(:);
        t=t(:);
        [P{i},~]=polyfit(t,s,PolyOrder);
        S=polyval(P{i},t);
        plot(t,S) 
    end
    AF=[];
else
if isequal(Pcs,'P,QRS,T,U')
    AF=ecgfidapprox(ecg,fs,hr,FAM); % finding approximated fiducial points
    AF2=AF; % making a copy of AF
    %%% changing the AF according to the ordered pieces    
    AF2(2,2)=AF(2,3);
    AF2(3,2)=AF(3,4);
    AF2(:,[3 4])=[]; np=size(AF2,2);
elseif isequal(Pcs,'P,QRS,T')
    AF=ecgfidapprox(ecg,fs,hr,FAM); % finding approximated fiducial points
    AF2=AF; % making a copy of AF
    %%% changing the AF according to the ordered pieces    
    AF2(2,2)=AF(2,3);
    AF2(3,2)=AF(3,4);
    AF2(:,[3 4 6])=[]; np=size(AF2,2);
elseif isequal(Pcs,'P,Q,R,S,T') 
    AF=ecgfidapprox(ecg,fs,hr,FAM); % finding approximated fiducial points
    AF2=AF; % making a copy of AF
    %%% changing the AF according to the ordered pieces    
    AF2(:,end)=[];
    np=size(AF2,2);
elseif isequal(Pcs,'P,Q,R,S,T,U')
    AF=ecgfidapprox(ecg,fs,hr,FAM); % finding approximated fiducial points
    AF2=AF; % making a copy of AF
    np=size(AF2,2);
elseif ~isequal(Pcs,'MovingWindow')
    error('wrong ecg pieces!')  
end    

%%% Ploting
figure
plot(ecg);
hold on
plot(AF2(:),ecg(AF2(:)),'r*');



%%% correcting PolyOrder vector if it is a scalar value
if size(PolyOrder(:),1)==1
    PolyOrder=PolyOrder.*ones(1,size(AF2,2));
elseif  size(PolyOrder(:),1)~=np
    error('wrong order of polynomials!')
end


%% Fitting the polynomials
%%% if the method is default
if nargin==4
    for i=1:np
        t=AF2(1,i):AF2(3,i);
        s=ecg(t);
        s=s(:);
        t=t(:);
        [P{i},~]=polyfit(t,s,PolyOrder(i));
        S=polyval(P{i},t);
        plot(t,S)   
    end

%%% if the method is fix R point    
elseif nargin==5 || isequal(PFM,'fixR')
    for i=1:np
        t=AF2(1,i):AF2(3,i);
        s=ecg(t);
        s=s(:);
        t=t(:);
            if AF(2,3)>t(1) && AF(2,3)<t(end)
            [P{i},~]=polyfix(t,s,PolyOrder(i),AF(2,3),ecg(AF(2,3)));
            else
            [P{i},~]=polyfit(t,s,PolyOrder(i));
            end
        S=polyval(P{i},t);
        plot(t,S)
    end

%%% if the method is fix three points Q, R and S    
elseif nargin==5 || isequal(PFM,'fixQRS')
switch Pcs
    case {'P,QRS,T,U','P,QRS,T'}
    for i=1:np
        t=AF2(1,i):AF2(3,i);
        s=ecg(t);
        s=s(:);
        t=t(:);
            if AF(2,3)>t(1) && AF(2,3)<t(end)
            [P{i},~]=polyfix(t,s,PolyOrder(i),AF(2,2:4),ecg(AF(2,2:4))');
            else
            [P{i},~]=polyfit(t,s,PolyOrder(i));
            end
        S=polyval(P{i},t);
        plot(t,S)
    end
    case {'P,Q,R,S,T','P,Q,R,S,T,U'}
        for i=1:np
        t=AF2(1,i):AF2(3,i);
        s=ecg(t);
        s=s(:);
        t=t(:);
           if i==2 || i==3 || i==4
            [P{i},~]=polyfix(t,s,PolyOrder(i),AF(2,i),ecg(AF(2,i))');
            else
            [P{i},~]=polyfit(t,s,PolyOrder(i));
            end
        S=polyval(P{i},t);
        plot(t,S)
        end     
end    
end
end 
end


