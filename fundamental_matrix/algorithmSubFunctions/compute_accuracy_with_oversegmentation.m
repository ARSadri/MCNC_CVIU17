function [percnt, prms] = compute_accuracy_with_oversegmentation(s, labels)
    out = zeros(max(s),max(labels));
    for scnt = 1 : max(s)
        for lcnt = 1 : max(labels)
            out(scnt,lcnt) = sum((s==scnt)+(labels==lcnt)==2);
        end
    end

    [~,lbl_new] = max(out);
    for idcnt = 1 : max(labels)
        labels(labels==idcnt) = -lbl_new(idcnt);
    end
    labels = -labels;  
    
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
