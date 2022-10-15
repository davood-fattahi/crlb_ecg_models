function x=GausVal(t,p)
p=p(:);NumGaus=length(p)/3;
g='@(p)';
for i=1:NumGaus
    g=([g '+(p(' num2str((i-1)*3+1) ')*exp(-((t-p(' num2str((i-1)*3+3) ')).^2)./(2*p(' num2str((i-1)*3+2) ')^2)))']);
end
g=eval(g);
x=g(p);

