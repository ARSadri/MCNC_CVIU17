function [G] = corr_H2G( H )
%find_corr finds the corelation between rows of a given H
% if H is MxN then G will be MxM

H = H-repmat(mean(H,2),1,size(H,2));
G = H*H'/size(H,2);
dA = eps+std(H,[],2);
G = G./(repmat(dA,1,size(G,2)).*repmat(dA',size(G,1),1));

end