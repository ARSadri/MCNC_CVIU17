function [H_T, R_T, sample_type, Lik] = MCMC_km_funda(X, model_type, label, numHypo)

%% args
% use labels to visualize where a sample is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_samples = numHypo;
n_FLKOS_itrs = 12;

Thresh_LJ = 2;
k_LJ = 12;
psize_LJ = 8;

Thresh_SJ = 2;
k_SJ = 20;
psize_SJ = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,n] = size(X);
[ fitfn, resfn] = getModelParam(model_type);

inds_old = randsample(n, psize_LJ, false);
psub = X(:,inds_old);
x_start = feval(fitfn,psub);
J_old = eps;
dist_old = feval(resfn,x_start,X);
SRes = sort(dist_old);
sigma_old = SRes(psize_LJ);
jmp_coin = 1;

H_T = zeros(n,n_samples);
R_T = H_T;
sample_type = zeros(1,n_samples);
Lik = zeros(1,n_samples);
totalr1 = 0;
totalr2 = 0;
totala1 = 0;
totala2 = 0;

smpls_cnt = 0;
while( (smpls_cnt < n_samples) && (totalr2+totalr1 < 6*n_samples) )
    if rand < jmp_coin
        %%%%%%%%%%%%%%%%%%%%%%%% X_start to X_PHI  %%%%%%%%%%%%%%%%%%%%%%%
%         u1 = hist(label(inds_old),1:3);
        x_start_aff = 1 - exp(- (0.2*dist_old/(2*sigma_old^2)).^4 );
        x_start_aff = x_start_aff/sum(x_start_aff);
        Next_Mode_k = max([psize_LJ randsample(n, 1, true, x_start_aff)]);
        [~,I] = sort(dist_old);
        inds = I(Next_Mode_k - psize_LJ+1: Next_Mode_k);
        psub = X(:,inds);
        x_PHI = feval(fitfn,psub);
        
        %%%%%%%%%%%%%%%%%%%%%%%% X_PHI to X_F  %%%%%%%%%%%%%%%%%%%%%%%
        iter_cnt = 0;
        while( iter_cnt < n_FLKOS_itrs )
            theta = feval(fitfn,psub);
            dist = feval(resfn,theta,X);
            dist(x_start_aff<0.1*max(x_start_aff))=inf;
            [~,I] = sort(dist);
            inds = I( k_LJ-psize_LJ+1:k_LJ );
            psub = X(:,inds);
            iter_cnt = iter_cnt+1;
        end
        dist_x_F = feval(resfn,theta,X);
        sigma_x_F = Fast_MSSE(dist_x_F.^0.5,Thresh_LJ, psize_LJ);
        if isinf(sigma_x_F)
            totalr1 = totalr1 +1;    
            continue;
        end
        aff_x_F = exp(-dist_x_F./(2*sigma_x_F^2));

        %%%%%%%%%%%%%%%%%%%%%%% X_F to Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        inds = I( k_LJ-psize_LJ+1:k_LJ );
        y_from_Q0 = feval(fitfn,X(:,inds));
        dist_y = feval(resfn,y_from_Q0,X);
        sigma_y = Fast_MSSE(dist_y.^0.5,Thresh_LJ, psize_LJ);
        if isinf(sigma_y)
            totalr1 = totalr1 +1;    
            continue;
        end
        SRes = sort(dist_y);
        J_y = exp(-SRes(k_LJ)/(2*sigma_y^2))/sigma_y;

        inds_new = inds;
%         u2=hist(label(inds),1:3);
        aff_y = exp(-dist_y./(2*sigma_x_F^2));
        dist_y_from_x_F = aff_x_F'*aff_y;
        
        %%%%%%%%%%%%%%%%%%%%%%% Y to Y_PHI %%%%%%%%%%%%%%%%%%%%
        PHI_F =  x_start\x_PHI;
        y_phi = y_from_Q0 / PHI_F;
        
        %%%%%%%%%%%%%%%%%%%%% Y_PHI to Y_F %%%%%%%%%%%%%%%%%%%%%%%
        theta = y_phi;
        iter_cnt = 0;
        while( iter_cnt < n_FLKOS_itrs )
            dist = feval(resfn,theta,X);
            [~,I] = sort(dist);
            inds = I( k_LJ-psize_LJ+1:k_LJ );
%             if iter_cnt==1
%                 u3=hist(label(inds),1:3);
%             end
            psub = X(:,inds);
            iter_cnt = iter_cnt+1;
            if( iter_cnt < n_FLKOS_itrs )
                theta = feval(fitfn,psub);
            end
        end
%         u4=hist(label(inds),1:3);
        dist_y_F = dist;
        
        %%%%%%%%%%%%%%%%%%%%%%% Q(X_Start|Y_F) %%%%%%%%%%%%%%%%%%%%%%
        aff_y_F = exp(-dist_y_F./(2*sigma_x_F^2));
        aff_x_start = exp(-dist_old./(2*sigma_x_F^2));
        dist_x_from_y_F = aff_y_F'*aff_x_start;
        
        %%%%%%%%%%%%%%%%%%%%%%% MetropolisHastings %%%%%%%%%%%%%%%%%%%%
        MH_crit = (J_y/J_old)*dist_x_from_y_F/dist_y_from_x_F;
        if isnan(MH_crit)
            totalr1 = totalr1 +1;    
            continue;
        end
        if  MH_crit > rand
%             disp(u1), disp(u2), disp(u3), disp(u4), disp(MH_crit)
            smpls_cnt = smpls_cnt + 1;
            J_old = J_y;
%             inds_old = inds_new;
            x_start = y_from_Q0;
            sample_type(smpls_cnt) = 1;
            Lik(smpls_cnt) = J_old;
            H_T(:,smpls_cnt) = exp(-dist_y/( 2*(sigma_y^2) ) );
            R_T(:,smpls_cnt) = dist_y.^0.5/sigma_y;
            totala1 = totala1 + 1;
            jmp_coin = mean(H_T(:,smpls_cnt));
            dist_old = dist_y;
            sigma_old = sigma_y;
        else
            totalr1 = totalr1 +1;    
        end
    else
        
        [~,I] = sort(dist);
        x_start_aff = exp(- (0.8*dist_old/(2*sigma_old^2)).^4 );        
        x_start_aff = x_start_aff/sum(x_start_aff);
        closeby_k = max([psize_LJ randsample(n, 1, true, x_start_aff)]);
        inds = I(closeby_k - psize_LJ+1: closeby_k);
        psub = X(:,inds);
        for itr = 1:n_FLKOS_itrs
            LS_theta = feval(fitfn,psub);
            dist = feval(resfn,LS_theta,X);
            [~,I] = sort(dist);
            inds_new = I(k_SJ-psize_SJ+1:k_SJ);
            psub = X(:,inds_new);
        end
        dist = feval(resfn,LS_theta,X);
        sigma = Fast_MSSE(dist.^0.5,Thresh_SJ, psize_SJ);
        if isinf(sigma)
            totalr2 = totalr2 +1;
            jmp_coin = 0.5;
            continue;
        end

        J_new = exp(-SRes(k_SJ)/(2*sigma_old^2))/sigma_old;
        M_crit = (J_new/J_old);

%         [u,~] = hist(label(inds_new),1:3); disp([u M_crit])

        if M_crit > rand
%             inds_old = inds_new;
            J_new = mean(exp(-dist/(2*sigma^2)))/sigma;
            smpls_cnt = smpls_cnt + 1;
            x_start = LS_theta;
            sample_type(smpls_cnt) = 2;
            H_T(:,smpls_cnt) = exp(-dist/( 2*(sigma^2) ) );
            R_T(:,smpls_cnt) = dist.^0.5/sigma;
            dist_old = dist;
            sigma_old = sigma;
            J_old = J_new;
            Lik(smpls_cnt) = J_old;
            jmp_coin = mean(H_T(:,smpls_cnt));
            totala2 = totala2 + 1;
        else
            jmp_coin = 0.5;
            totalr2 = totalr2 +1;
        end

    end
end
disp([totala1 totala2; totalr1 totalr2;...
    floor(100*(totala1+totala2)/(totala1+totala2+totalr1+totalr2)) ...
    floor(100*(totalr1+totalr2)/(totala1+totala2+totalr1+totalr2))  ])
close all    
end