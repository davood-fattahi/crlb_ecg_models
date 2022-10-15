
function PP=deviate(P,pdev)
for i=1:size(P,1)
    ppp=flip(P{i,1});
    Pdev=zeros(size(P{i,1}));
    for m=1%size(P{l,1},2)
        Pdev(m)=(pdev)^(m).*randn(1);
    end 
    PP{i,1}=P{i,1}+flip(Pdev);
end
end

