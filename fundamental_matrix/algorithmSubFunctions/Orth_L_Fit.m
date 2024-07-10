function u = Orth_L_Fit( x,y )
%Orthogonal distance based line fitting
% inputs: x and y are row vectors of data points
% outputs: u has two elements, slope and intercept

n = length(x);
uy = mean(y);
ux = mean(x);
B = 0.5*((sum(y.^2) - n*uy^2)-(sum(x.^2) - n*ux^2))/(n*ux*uy-sum(x.*y));
b1 = -B + (B^2 + 1)^0.5;
a1 = uy - b1*ux;
u1 = [b1 a1];
res1 = [-u1(1) 1 -u1(2)]*[x;y;1+0*x]/(norm(u1));
b2 = -B - (B^2 + 1)^0.5;
a2 = uy - b2*ux;
u2 = [b2 a2];
res2 = [-u2(1) 1 -u2(2)]*[x;y;1+0*x]/(norm(u2));
if mean(abs(res1)) < mean(abs(res2))
    u = u1;
else
    u = u2;
end

end

