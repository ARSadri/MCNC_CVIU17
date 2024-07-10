function [dist] = motion_res(P, X)
[F_2,N] = size(X);
uX = P(:,end);
X = X - repmat(uX,1,N);
P = P(:,1:end-1);
dist = sum((P*P'*X - X).^2)'/N^2;
TST = ones(F_2,1);
distTST = sum((P*P'*TST - TST).^2);
dist = dist/distTST;

end