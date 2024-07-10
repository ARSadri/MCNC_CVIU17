function res = line1D_res(P,X)

x = X(1,:)';
y = X(2,:)';
n = length(x);

res = (y-[x]*P).^2;
