function [ lower, matdiag, upper ] = lowerdiagupper( A )
% lowerdiagupper accepts a square matrix and gives back
% three parts of it in form of vectors
% lower tri part       diagonal part      upper tri part
upper = zeros((length(A)*length(A) - length(A))/2,1);
lower = upper;
matdiag = zeros(length(A),1);
lcnt = 0;
ucnt = 0;
dcnt = 0;
for k = 1 : length(A)
    for j = 1 : length(A)
        if(k < j)
            ucnt = ucnt + 1;
            upper(ucnt) = A(k,j);
        elseif(k > j)
            lcnt = lcnt + 1;
            lower(lcnt) = A(k,j);
        else
            dcnt = dcnt + 1;
            matdiag(dcnt) = A(k,j);
        end
    end
end

end

