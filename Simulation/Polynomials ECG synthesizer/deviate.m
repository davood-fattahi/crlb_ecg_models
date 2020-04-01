
function PP=deviate(P,pdev)
for l=1:size(P,1)
    ppp=flip(P{l,1});
    Pdev=zeros(size(P{l,1}));
    for m=2%size(P{l,1},2)
        Pdev(m)=(pdev)^(m).*randn(1);
    end 
    PP{l,1}=P{l,1}+flip(Pdev);
end
end

