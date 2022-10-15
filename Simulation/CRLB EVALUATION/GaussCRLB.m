function [CRLB, FI]=GaussCRLB(a, b, Fs, varn )

Iaa=Fs*b*sqrt(pi);
Iab=Fs*a*sqrt(pi)./2;
Iac=0; Ibc=0;
Ibb=(Fs*a^2*3*sqrt(pi))./(b*4);
Icc=(Fs*a^2*sqrt(pi))./(b*2);
FI=[Iaa Iab Iac; Iab Ibb Iac; Iac Ibc Icc]./varn;
CRLB=varn.*inv([Iaa Iab Iac; Iab Ibb Iac; Iac Ibc Icc]);


