function [ s,theta_m, Tx, Ty ] = Affine_motion_mdoel_RANSAC_Hopkins(clus)
%Affine_motion_mdoel gives out the scale, Rotation and the start and end of
%translation
%clus must be 3-tensor, number_of_points - x,y - number_of_frames
%intv_strt begingn frame
%intv_end  last frame

    k_order = min([ 10 max([ceil(size(clus,1)*0.03) 2])]);
    if k_order>10
        k_order = 10;
    end
    if k_order < 2
        k_order=2;
    end
    no_samples = min([200 max([2 1+ceil(size(clus,1)*0.4)])]);
        k_order=min([10 no_samples]);

    inds = zeros(k_order,no_samples);
    for cnt = 1 : no_samples
        inds(:,cnt) = randsample(size(clus,1),k_order);
    end
    struct_error = 10^10*ones(1,no_samples);
    for smplcnt = 1 : no_samples
        CC = clus(inds(:,smplcnt), :, 1)';
        CC1 = clus(inds(:,smplcnt), :, end)';
        rlt0 = mean(CC, 2); 
        rlt10 = mean(CC1, 2); 
        rlt = CC - repmat(rlt0, [1, size(CC, 2)]); 
        rlt1 = CC1 - repmat(rlt10, [1, size(CC1, 2)]); 
        s = sqrt(max(rlt1(:).^2) / max(rlt(:).^2)); 
        theta = zeros(1,size(rlt, 2));
        for i=1:size(rlt, 2)
            a = rlt(:, i);
            b = rlt1(:, i);
            theta(i) = acos(dot(a,b)/(norm(a)*norm(b)));
        end
        theta(isnan(theta))=0;
        theta_m = mean(theta);
        Ah = [0 -1
              1  0];
        R = eye(2) + sin(theta_m)*Ah + (1-cos(theta_m))*Ah^2;
        CC = clus(:, :, 1)';
        CC1 = clus(:, :, end)';

        CC_rlt = CC - repmat(rlt0, [1, size(CC, 2)]); 
        CC_rlt1 = CC1 - repmat(rlt10, [1, size(CC, 2)]); 
        u1 = s*R*CC_rlt;
        u2 = CC_rlt1;
        u1n = sum(u1.^2).^0.5;
        u2n = sum(u2.^2).^0.5;
        uns = sort([u1n;u2n]);
        cos_theta = sum(u1.*u2)./(u1n.*u2n);
        motmodediff = 0.5*(1-cos_theta.*uns(1,:)./uns(2,:));
        udiffn = sum((u2-u1).^2).^0.5;
        dist = motmodediff.*udiffn.*(uns(2,:).^2);
        dist(isnan(dist))=10^10;
        the_dist = dist(dist<10^9);
        if ~isempty(the_dist)
            struct_error(smplcnt) = mean(the_dist);
        end
    end
    [~,bestind] = min(struct_error);
    bestind = bestind(1);
    
    CC = clus(inds(:,bestind), :, 1)';
    CC1 = clus(inds(:,bestind), :, end)';
    rlt0 = mean(CC, 2); 
    rlt10 = mean(CC1, 2); 
    rlt = CC - repmat(rlt0, [1, size(CC, 2)]); 
    rlt1 = CC1 - repmat(rlt10, [1, size(CC1, 2)]); 
    
    s = sqrt(sum(rlt1(:).^2) / sum(rlt(:).^2)); 
    theta = zeros(1,size(rlt, 2));
    for i=1:size(rlt, 2)
        a = rlt(:, i);
        b = rlt1(:, i);
        theta(i) = acos(dot(a,b)/(norm(a)*norm(b)));
    end
    theta_m = mean(theta);
    T=rlt10-rlt0;
    Tx = T(1);
    Ty = T(2);
end