function [hstd,INLvec, Len] = MSSE( fiterrs, Lambda, m_p, k_m)
    Lambda_sq = Lambda^2;
    INLvec = 0*fiterrs;    
    hstd = inf;
    Len = 0;
    numpts = length(fiterrs);
    %structure should be bigger than 1% and smaller than 99% of data points
    struct_min_size = max([m_p+1 k_m ceil(numpts/100)]);
    struct_max_size = numpts - ceil(numpts/100);
    if numpts < struct_min_size
        return
    end
    
    [R2, R2ind]= sort(fiterrs);
    msse_crit = ((1 : numpts)'-m_p).*R2./cumsum(R2);
    msse_crit(1:struct_min_size)=0;
    msse_crit(msse_crit<Lambda_sq)=0;
    [~, ptcnt] = max(msse_crit>0);
    
    if (ptcnt > struct_min_size)  && (ptcnt < struct_max_size)
        INLvec(R2ind(1:ptcnt-1))=1;
        hstd = (sum(R2(1:ptcnt))/(ptcnt-1-m_p))^0.5;
        Len = sum(INLvec==1);
    end
end