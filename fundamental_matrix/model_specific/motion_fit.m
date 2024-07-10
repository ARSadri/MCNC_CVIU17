function [P] = motion_fit(Xs)

uX = mean(Xs,2);
X1 = Xs - repmat(uX,1,size(Xs,2));
[U,sV,~] = svd(X1,'econ');

sV = sV/size(Xs,2)^0.5;

sV = diag(sV);
sV = sV./[1;cumsum(sV(1:end-1))];
[~, d] = max(sV<0.01);
d = min(d,4);
if isempty(d)
    d = 4;
end

Y = U(:,1:d);
P = [Y uX];
% mean(sum((Y*Y'*X1-X1).^2).^0.5)

% matrank = modelselection(diag(s),100e-6);
% d = max(2,matrank);
% d = min(4,d);

% disp([s d-1])
% Y = U(:,1:d-1);
% P = [Y uX];

% d=4;
% s = diag(s);
% ind = find(s<1.0,1,'first');
% if ind >4
%     ind = 4;
% end
% d=4;
% Y = U(:,1:r);
% H = eye(size(X,1)) - Y*Y';
% P = [Y*Y' uX];%Y*((Y'*Y)\Y');
% P = H(:);

% Y = U(:,1:d);
% P = [Y uX];

end