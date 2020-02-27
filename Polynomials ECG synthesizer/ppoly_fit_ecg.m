function [P, AprxSegs] =ppoly_fit_ecg(ecg,fs,PolyOrder,np, varargin)


hr=fs./size(ecg(:),1);

AprxSegs=EcgSegmentorApprox(ecg,fs,hr,np);
AprxSegs=AprxSegs(:);

figure
plot(ecg);
hold on
plot(AprxSegs(:),ecg(AprxSegs(:)),'r*');
    
if nargin==4
    for i=1:np
        t=AprxSegs(2*(i-1)+1):AprxSegs(2*i);
        s=ecg(t);
        s=s(:);
        t=t(:);
        if size(PolyOrder(:),1)==np
        [P,~]=polyfit(t,s,PolyOrder(i));
        elseif size(PolyOrder(:),1)==1
        [P,~]=polyfit(t,s,PolyOrder);
        else
        error('wrong order of polynomials')
        end
        S=polyval(P,t);
        plot(t,S)   
    end
  

elseif nargin==5 || varargin{1}=='fixR'
    [Ra, Ri]=max(ecg);
    tt=AprxSegs(5):AprxSegs(6);
    B=(ecg(tt)-Ra); B=B(:)';
    A=((tt-Ri).^2); A=A(:)';
    b=B*A'*((A*A')^(-1));
    P=[b -2*b*Ri b*Ri^2+Ra];
    S=polyval(P,tt);
    plot(tt,S)  
    
end


