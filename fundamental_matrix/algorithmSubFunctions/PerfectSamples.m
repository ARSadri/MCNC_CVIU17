function [H_T, R_T, Likli] = PerfectSamples(X, model_type, G, n_samples, labels)

[fitfn, resfn,~,psize,] = getModelParam(model_type);

for lbcnt = 1 : max(labels)
    inds = labels==lbcnt;
    psub = X(:,inds);
    theta = feval(fitfn,psub);
    dist = feval(resfn,theta,X);
    [~,I] = sort(dist);
    sigma = Fast_MSSE(dist.^0.5,2.15, 6);        
    aff = exp(-dist./(2*sigma^2));
    J_y = aff(I(12))/sigma;
    Likli(lbcnt) = J_y;
    H_T(:,lbcnt) = aff;
    R_T(:,lbcnt) = dist.^0.5/sigma;
end

end