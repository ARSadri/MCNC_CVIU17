function [CutPoint,stdev] = applyMSSE(r2,k,T)

if(k>= length(r2))
    CutPoint = length(r2);
else
    sig1 = cumsum(r2)./[1;1;(1:length(r2)-2)'];

    sig = sig1(k:end-1);
    r = r2(k+1:end);

    mask = (r>sig*T*T);

    fInd = find(mask>0,1,'first') ;

    if(isempty(fInd))

%             fInd = 0;
%     elseif(fInd > 0.9*length(sig))
       fInd = length(r2)-k;    
    end

    fInd = fInd+k-1;

    CutPoint = fInd;
    stdev = sig1(fInd)^.5;


end
