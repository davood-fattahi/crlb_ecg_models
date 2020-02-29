function [P, WS, PO, NN]=polyfitAMW(t,s,tr,maxi,varargin)

% 
% [P, S]=Apolyfit(t,s,eps,'PolyOrder',po,wsi}
% [P, S]=Apolyfit(t,s,eps,'WindowSize',ws,poi}


t=t(:);
s=s(:);
s=normalize(s,1,'scale');
    

n=1:size(s(:),1);
switch varargin{1}
    case 'PolyOrder'
        po=varargin{2};
        if isempty(varargin{3})
        wsh=size(s(:),1);
        elseif size(varargin{3}(:),1)==1
        wsh=varargin{3};
        wsl=2*po+1;
        elseif size(varargin{3}(:),1)==2
        wsh=varargin{3}(2);
        wsl=varargin{3}(1);
        end
        
        for j=0:po
        vc(j+1)=nchoosek(po,j);
        end
        vc=sum(vc.^2);       
        
        nn=n(1);
        WS=[];
        NN=[];
        PO=[];
        while nn<n(end)    
            i=1;
            if nn+wsh-1<n(end)
                ws=wsh;
            else
                ws=n(end)-nn+1;
            end
            if ws>wsl
            ss=s(nn:nn+ws-1);
            D=diff(ss,po);
            ee=var(D)./vc;
            while ee>tr && i<maxi && ws>wsl
                ws=ws-5;
                D(end-5+1:end)=[];
                ee=var(D)./vc;
                i=i+1; 
            end
            end
            WS=[WS ws];
            NN=[NN nn];
            nn=nn+ws;
            PO=[PO po];
        end
        
        
        
    case 'WindowSize'
        ws=varargin{2};
        if isempty(varargin{3})
        poi=6;
        else
        poi=varargin{3};
        end
        nn=n(1);
        WS=[];
        NN=[];
        PO=[]; 
        while nn<n(end)    
            i=1;
            if nn+wsh-1<n(end)
                ws=wsh;
            else
                ws=n(end)-nn+1;
            end
            ss=s(nn:nn+ws-1);
            D=diff(ss,po);
            ee=abs(mean(D.^2,'all')-max(D.^2,[],'all'));
            while ee>tr && i<maxi
                po=po+1;
                D=diff(ss,po);
                ee=abs(mean(D.^2,'all')-max(D.^2,[],'all'));
                i=i+1;
            end
            WS=[WS ws];
            NN=[NN nn];
            nn=nn+ws;
            PO=[PO po];
        end 
        
        
end

figure
plot(t,s)
for i=1:size(NN(:),1)
    nn=NN(i):NN(i)+WS(i)-1;
    tt=t(nn); tt=tt(:);
%% LS:
    A=ones(WS(i),PO(i)+1);
    for j=PO(i):-1:1
        A(:,j)=A(:,j+1).*tt;
    end
    [Q,R] = qr(A,0);
    P{i} = R\(Q'*s(nn));

%%% Or:
%     P{i}=lsqr(A,s(nn));

%%% Or:
%     [P{i},~]=polyfit(tt,s(nn),PO(i));


%%%% Or
%     P{i}=((A'*A)^-1)*A'*s(nn);
    
    S=polyval(P{i},tt);
    hold on
    plot(t(nn),S)
end
end        
        