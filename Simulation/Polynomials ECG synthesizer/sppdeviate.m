function varargout=sppdeviate(pp,beta,gamma,amp)

% dpp=sppdeviate(pp,beta,gamma,amp)


%%

coefs=pp.coefs;
breaks=pp.breaks;

% de-localization of the coefs
for i=1:size(coefs,1)
    coefs(i,:) = polydeviate(coefs(i,:),1,-breaks(i),1);
end

% deviate coefs
[coefs,breaks] = polydeviate(coefs,beta,gamma,amp,breaks);

% re-localize the coefs
for i=1:size(coefs,1)
    coefs(i,:) = polydeviate(coefs(i,:),1,breaks(i),1);
end

% deliver the outputs
pp.coefs=coefs;
pp.breaks=breaks;
varargout{1}=pp;
end
