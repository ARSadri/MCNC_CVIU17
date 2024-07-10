function [ filven ] = FastAvgFilter( invec, filsize )

    filsize = max([2 2*floor(filsize/2)]);
    filven = filsize*max(invec)*ones(1,length(invec));
    filven(filsize) = sum(invec(1:filsize));
    for ptcnt = filsize+1:length(invec)
        filven(ptcnt) = filven(ptcnt-1) + invec(ptcnt) - invec(ptcnt-filsize);
    end
    filven = [filven(1:filsize/2+1) filven(filsize+1 : end) filven(filsize/2+1:filsize-1)] ;

end

