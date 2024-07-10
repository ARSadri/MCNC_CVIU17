function [ R,x,theta ] = RotAtoB( a,b )
%RotAtoB Gives the rotation matrix that rotates the vector
% b to a using Rodrigues' rotation formula, R = I + (Sin t)A + (1-cos t)A^2
    %x = [a(2)*b(3) - b(2)*a(3);a(3)*b(1) - b(3)*a(1);a(1)*b(2) - b(1)*a(2)];
    x = cross(a,b);
    x = x/norm(x);
    theta = acos(dot(a,b)/(norm(a)*norm(b)));
    A = [0    -x(3)  x(2)
         x(3)   0   -x(1)
        -x(2)  x(1)   0];
    R = eye(3) + sin(theta)*A + (1-cos(theta))*A^2;
end

