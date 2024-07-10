function [ hestd, inc_vec ] = Fast_AVG_MSSE( fit_errs, Lambda, MSSE_max_num_cl, m_p, filsize )
%Fast_AVG_MSSE assumes that input residuals have a combination of 
%   uniform noise and guassian.
%   we first run an avg filter on data points and then run MSSE
%   hence the new parametr, filsize for avg filter.
% default : filsize>2
%           2<Lambda<3
%           MSSE_max_num_cl = 100   smallest structre has 1% of data points.
%           m_p number of model parameters
    fit_errs = fit_errs.^.5;

    filsize = 2*floor(filsize/2);
    Lambda_sq = Lambda^2;
    inc_vec = 0*fit_errs;    
    hestd = 0;    

    [x, xind]= sort(fit_errs);
    struct_minsize = max([ 20 ceil(length(fit_errs)/MSSE_max_num_cl)]);
    struct_maxsize = length(fit_errs) - ceil(length(fit_errs)/MSSE_max_num_cl);
    if length(x)<struct_minsize
        return
    end

    x = x.^2;
    x_f = x;
    S2 = sum(x_f(1:struct_minsize));
    for ptcnt = struct_minsize + 1 : length(x) - filsize/2
        x_f(ptcnt) = mean(x(ptcnt - filsize/2:ptcnt + filsize/2));
        MSCr(ptcnt) = x_f(ptcnt)*((ptcnt-1-m_p)/S2);
        if MSCr(ptcnt)>=Lambda_sq
            break;
        end
        S2 = S2 + x_f(ptcnt);
    end
    if ptcnt < struct_maxsize
        hestd = std(x(1:ptcnt).^0.5);
        inc_vec(xind(1:ptcnt))=1;
    end
    
    inc_vec = find(inc_vec==1);
end