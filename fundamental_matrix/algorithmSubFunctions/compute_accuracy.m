function [percnt, prms] = compute_accuracy(s, labels)
    P = perms(1:max(s)); 
    mn = zeros(1, size(P, 1)); 
    for i=1:size(P, 1)
        prms = P(i, :); 
        mn(i) = sum(abs(sign(s - prms(labels)))); 
    end
    [min_mn, id] = min(mn);
    percnt = 100*(numel(s) - min_mn)/numel(s);
    prms = P(id, :);
end
