function x=GausVal2(t,p)
% demo code for evaluation the model of (Gaussian+Line)
NumGaus=length(p)/5;
g='@(p)';
for i=1:NumGaus
    g=([g '+(p('  num2str((i-1)*5+1) ')*(t-'  num2str((i-1)*5+2) ')+p(' num2str((i-1)*5+3) ').*exp(-((t-p(' num2str((i-1)*5+5) ')).^2)./(2*p(' num2str((i-1)*5+4) ')^2)))']);
end
g=eval(g);
x=g(p);

