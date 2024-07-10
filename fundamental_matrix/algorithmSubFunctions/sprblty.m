function [ out ] = sprblty( fit_errs, Lambda1, Lambda2 )

    x= sort(fit_errs);
    out = zeros(length(x));
    for ptcnt_f = 1 : length(x) - 1
        out(ptcnt_f,1:ptcnt_f) = x(ptcnt_f)/std(x(1:ptcnt_f));
        for ptcnt = ptcnt_f + 1 : length(x)
            out(ptcnt_f,ptcnt) = mean(x(ptcnt_f+1:ptcnt))/(Lambda1*std(x(1:ptcnt_f))+Lambda2*std(x(ptcnt_f+1:ptcnt)));
        end
    end      

end