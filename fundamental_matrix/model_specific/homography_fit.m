function P = homography_fit(X)

x1 = X(1:3,:);
x2 = X(4:6,:);

H = vgg_H_from_x_lin(x1,x2);
%[H,~] = vgg_H_from_x_nonlin(H,x1,x2);
if isempty(H)
    H=zeros(3);
end
P = H(:);

end


