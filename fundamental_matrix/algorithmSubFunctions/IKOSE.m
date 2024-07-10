function hstd = IKOSE( fiterrs, Lambda)
    hstd = inf;
    numpts = length(fiterrs);
    %structure should be bigger than 1% and smaller than 99% of data points
    struct_min_size = max([12 ceil(numpts/100)]);
    struct_max_size = numpts - ceil(numpts/100);
    if numpts < struct_min_size
        return
    end
    K = struct_min_size;
    err_srtd = sort(fiterrs);
    sigma = max(fiterrs(fiterrs<inf));
    for itr = 1 : 10
        n = sum(err_srtd<Lambda*sigma);
        Kapa = K/n;
        sigma = err_srtd(K)/norminv((1+Kapa)/2,0,1);
    end
    if (n > struct_min_size)  && (n < struct_max_size)
        hstd = err_srtd(n)/Lambda;
    end
end