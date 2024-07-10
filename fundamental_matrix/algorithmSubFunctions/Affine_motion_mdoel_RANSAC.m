function [ s,R,rlt0,rlt10 ] = Affine_motion_mdoel_RANSAC(clus)
% This function rendomly uniformly selects some k_order samples from the structure
% and generates a least square model and finds fitting error of all data
% points of the structure with the model
% then choses the model with least median fitting error


    k_order = min([ 10 max([ceil(size(clus,1)*0.03) 2])]);
    if k_order>10
        k_order = 10;
    end
    if k_order < 2
        k_order=2;
    end
    no_samples = min([500 max([2 1+ceil(size(clus,1)*0.1)])]);
        k_order=min([10 no_samples]);

    inds = zeros(k_order,no_samples);
    for cnt = 1 : no_samples
        inds(:,cnt) = randsample(size(clus,1),k_order);
    end
    inds = unique(sort(inds,1)','rows')';
    inds(:,sum(diff(inds)==0,1)>0)=[];  %samples are not degenerate
    no_samples = size(inds,2);
    
    struct_error = 10^10*ones(1,no_samples);
    for smplcnt = 1 : no_samples
        
        [ s,R,rlt0,rlt10 ] = Affine_motion_mdoel(clus(inds(:,smplcnt), :, [1 end]));
        
        CC = clus(:, :, 1)';
        CC1 = clus(:, :, end)';

        CC_rlt = CC - repmat(rlt0, [1, size(CC, 2)]); 
        CC_rlt1 = CC1 - repmat(rlt10, [1, size(CC, 2)]); 
        u1 = s*R*CC_rlt;
        u2 = CC_rlt1;
        
        dist = sum((u2-u1).^2).*sum(CC_rlt1.^2).^0.5;
        dist(isnan(dist))=10^10;
        the_dist = dist(dist<10^9);
        if ~isempty(the_dist)
            struct_error(smplcnt) = median(the_dist);
        end
    end
    [~,bestind] = min(struct_error);
    bestind = bestind(1);
    
    [ s,R,rlt0,rlt10 ] = Affine_motion_mdoel(clus(inds(:,bestind), :, [1 end]));

end