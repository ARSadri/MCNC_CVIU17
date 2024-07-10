function [ s,R,rlt0,rlt10 ] = Affine_motion_mdoel(clus)
%Affine_motion_mdoel calculates affine transform parameters for input
%cluster
%   the affine motion model includes a 3-translation, 1-rotation and
%   1-scale all together 5 parametrs for N matching data points seen
%   through two frames.
% Inputs: clus : input cluster of points should be a N by 2 by 2 matrix
% N being the number of data points, then the first 2 is for x and y
% coordinates of every data point and the second 2 will be the 1st frame
% and the 2nd frame
% Outputa: s:scale
%          R: 2x2 Rotation matrix using lie algebra , 2nd order taylor
%          rlt0 and rlt10 : means of corrdinates of points in two frames
%          Translation vector is rlt10-rlt0

%Affine motion model has two scales and a sheer but we ignore those

    s_sigma = 0.015;
    theta_m_sigma = .005;


    CC = clus(:, :, 1)';
    CC1 = clus(:, :, end)';
    
    % Translation
    rlt0 = mean(CC, 2); 
    rlt10 = mean(CC1, 2); 
    rlt = CC - repmat(rlt0, [1, size(CC, 2)]); 
    rlt1 = CC1 - repmat(rlt10, [1, size(CC1, 2)]); 

    % Scale
    s = mean(sum(rlt1.^2).^0.5./sum(rlt.^2).^0.5); 
    s = 1+(s-1).*(1-exp(-(s-1).^2/(2*s_sigma^2)));
    
    % Rotation using Rodrigues and average of rotation angle of all points
    % around camera principle  - Z - axis
    theta = zeros(1,size(rlt, 2));
%     for i=1:size(rlt, 2)
%         a = [rlt(:, i);0];
%         b = [rlt1(:, i);0];
% %         expvec = cross(a,b)/(norm(a)*norm(b));
% %         expvec = [0;0;a(1)*b(2)-a(2)*b(1)]/(norm(a)*norm(b));
%         expvec = [0;0;a(1)*b(2)-a(2)*b(1)]/((a(1)^2+a(2)^2)^0.5*(b(1)^2+b(2)^2)^0.5);
%         theta(i) = asin(expvec(3));
%     end
    
    for i=1:size(rlt, 2)
        a = rlt(:, i);
        b = rlt1(:, i);
        theta(i) = asin((a(1)*b(2)-a(2)*b(1))/((a(1)^2+a(2)^2)^0.5*(b(1)^2+b(2)^2)^0.5));
    end
    
    theta(isnan(theta))=0;
    theta_m = mean(theta);
    theta_m = theta_m.*(1-exp(-theta_m.^2/(2*theta_m_sigma^2)));
    Ah = [0 -1
          1  0];
    R = eye(2) + sin(theta_m)*Ah + (1-cos(theta_m))*Ah^2;
end