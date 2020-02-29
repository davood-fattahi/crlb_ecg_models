function k=ecgknots(ecg,depth,tr)



k=nan(2^depth,2);
k(1,1)=1;
k(1,2)=size(ecg(:),1);

plot(ecg)
hold on
for j=1:depth
for i=1:2^(j-1)
    s=ecg(k(i,1):k(i,2));
    line=s(k(i,1))+(1:size(s(:),1)).*((s(k(i,2))-s(k(i,1)))/size(s(:),1));
    hold on
    plot(line)
    [M,I]=max(s-line);
    k(i+1,2)=k(i,2);
    k(i,2)=I;
    k(i+1,1)=I;
end
end