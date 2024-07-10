function [H_T, R_T, Likli] = MCNC_funda(X, model_type, G, n_samples, labels)
[~,sind]=sort(labels);

n_FLKOS_itrs = 8;
k_m = 8;
MSSEThresh = 2.35;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[F_2, N] = size(X);

if isnan(G)
    G = zeros(N);
end

a = 1:N; a = a'; pos = 0*a; pos(sind)=a; clear a;

H_T = zeros(N,n_samples);
R_T = H_T;
Likli = zeros(1,n_samples);

[fitfn, resfn,~,psize,] = getModelParam(model_type);
psize = 2+psize;

psize = 8;

%Intialization by a random state
sigma_old = inf;
while isinf(sigma_old)
    inds = randsample(N, psize, false);
%     disp(labels(inds))
    psub = X(:,inds);
    iter_cnt = 0;
    while( iter_cnt < n_FLKOS_itrs )
        theta = feval(fitfn,psub);
        dist_old = feval(resfn,theta,X);
        [~,I_old] = sort(dist_old);
        inds = I_old( k_m-psize+1:k_m);
%         disp(labels(inds))
        psub = X(:,inds);
        iter_cnt = iter_cnt+1;
    end
    sigma_old = MSSE(dist_old, MSSEThresh, psize, k_m);
end
J_old = eps;
jmp_coin = 1;
postBasedK_old = k_m;
psize_tmp_old = psize;
x_start = theta;

Num_reject_1 = 0;
Num_reject_2 = 0;
Num_reject_3 = 0;
Num_reject_4 = 0;
Num_reject_5 = 0;
Num_reject_6 = 0;
Num_reject_7 = 0;

% figure
% psub = X(:,labels<=1);
% theta = feval(fitfn,psub);
% dist = feval(resfn,theta,X);
% sigma = MSSE(dist, MSSEThresh, psize, k_m);
% subplot(2,3,1), plot(dist(sind).^0.5,'*'), title(sigma)
% subplot(2,3,4), plot(dist(sind).^0.5/sigma,'*')
% psub = X(:,labels==2);
% theta = feval(fitfn,psub);
% dist = feval(resfn,theta,X);
% sigma = MSSE(dist, MSSEThresh, psize, k_m);
% subplot(2,3,2), plot(dist(sind).^0.5,'*'), title(sigma)
% subplot(2,3,5), plot(dist(sind).^0.5/sigma,'*')
% psub = X(:,labels==3);
% theta = feval(fitfn,psub);
% dist = feval(resfn,theta,X);
% sigma = MSSE(dist, MSSEThresh, psize, k_m);
% subplot(2,3,3), plot(dist(sind).^0.5,'*'), title(sigma)
% subplot(2,3,6), plot(dist(sind).^0.5/sigma,'*')

% inds_rec = [];

MAXIMUM_EXPLORATION = 10*n_samples;

smpls_cnt = 0;
Total_jumps = 0;
while( (smpls_cnt < n_samples) && (Total_jumps-smpls_cnt < MAXIMUM_EXPLORATION) )
    
    Total_jumps = Total_jumps + 1;

    if rand < jmp_coin
        %%%%%%%%%%%%%%%%%%%%%%%% X_start to X_PHI  %%%%%%%%%%%%%%%%%%%%%%%
        x_start_aff = exp(-dist_old/(2*sigma_old^2));
        [~,I] = sort(dist_old);
        pts_weights = 1 - x_start_aff(I);
        pts_weights(pts_weights<0.95)=0;
        
%         Gp = G(I,I);
%         prob = zeros(N,1);
%         Gmax = max(Gp)';
%         for smplcnt = 1 : N - psize + 1
%             GSub = Gp(smplcnt : smplcnt + psize - 1,smplcnt : smplcnt + psize - 1);
%             prob(smplcnt) = mean(GSub(:))/max(Gmax(smplcnt : smplcnt + psize - 1));
%         end
%         prob(isnan(prob))=0.5;
%         if max(prob)==0
%             prob = 1 + prob;
%         end
%         pts_weights = pts_weights.*prob;
        
        Next_Mode_k = max([psize randsample(N, 1, true, eps+pts_weights)]);
        
        inds = I(Next_Mode_k - psize+1: Next_Mode_k);
        psub = X(:,inds);        
        x_PHI = feval(fitfn,psub);
        dist_x_PHI = feval(resfn,x_PHI,X);
        [~,I] = sort(dist_x_PHI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        GI = G(I,I);
        Pts2HypAvgp = zeros(1,N);
        Pts2HypAvgp(1:k_m)=inf;
        subG1sum = sum(sum(GI(1:k_m,1:k_m)));
        subG2sum = sum(sum(GI(k_m+1:N,k_m+1:N)));
        cutG = sum(sum(GI(1:k_m,k_m+1:N))) + sum(sum(GI(k_m+1:N,1:k_m)));
        for kcnt = k_m + 1 : N-1
            vec = GI(:,kcnt);
            to_add_to_1 = 2*sum(vec(1:kcnt)) - vec(kcnt);
            to_take_from_2 = 2*sum(vec(kcnt:N)) - vec(kcnt);
            subG1sum = subG1sum + to_add_to_1;
            subG2sum = subG2sum - to_take_from_2;
            cutG = cutG - to_add_to_1 + to_take_from_2;
            cutCost = cutG/min([subG1sum subG2sum]);
            if ( (cutCost>1.5*min(Pts2HypAvgp(1:kcnt-1))) && (cutCost<0.5*Pts2HypAvgp(k_m+1)) ) || min([subG1sum subG2sum])==0
                break;
            end
            Pts2HypAvgp(kcnt) = cutCost;
        end
        Pts2HypAvgp(1:k_m)=0;
        Pts2HypAvgp(isnan(Pts2HypAvgp))=0;
        Pts2HypAvgp(Pts2HypAvgp<0)=0;

        postBasedK = 12;
        if max(Pts2HypAvgp)>0
            postBasedK = randsample(N,1,true,Pts2HypAvgp);
        end
        postBasedK = min([postBasedK 40]);
        psize_tmp = max([psize floor(postBasedK/5)]);
        
        for  iter_cnt = 1 : n_FLKOS_itrs
            inds = I( postBasedK-psize_tmp+1:postBasedK );
            psub = X(:,inds);
            theta = feval(fitfn,psub);
            dist = feval(resfn,theta,X);
            [~,I] = sort(dist);
        end
        dist_x_F = dist;
%         clf, plot(dist_x_PHI(sind),'*'), hold on, plot(dist_x_F(sind),'*')
        sigma_x_F = MSSE(dist_x_F, MSSEThresh, psize_tmp, postBasedK);
        if isinf(sigma_x_F)
            Num_reject_1 = Num_reject_1 +1;
            continue;
        end
        aff_x_F = exp(-dist_x_F./(2*sigma_x_F^2));
        
        %%%%%%%%%%%%%%%%%%%%%%% X_F to Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        inds = I( postBasedK-psize_tmp+1:postBasedK );
        y_from_Q0 = feval(fitfn,X(:,inds));
        dist_y = feval(resfn,y_from_Q0,X);
%         subplot(3,1,3), plot(dist_y(sind).^0.5,'*')
        [~,I]= sort(dist_y);
        [sigma_y, ~, Len] = MSSE(dist_y,MSSEThresh, psize_tmp, postBasedK);
        if isinf(sigma_y)
            Num_reject_2 = Num_reject_2 +1;    
            continue;
        end
%         maxdist = MaxSqRes^0.5;
        J_y = exp(-mean(dist_y(I(1:Len)))/(2*(Len^0.5*sigma_y)^2))/(Len^0.5*sigma_y);
%         J_y = exp(-dist_y(I(postBasedK))/(2*sigma_y^2))/sigma_y;
%         J_y = (Len/(sigma_y*maxdist)) * ((normcdf(maxdist/sigma_y) - normcdf(-maxdist/sigma_y))^(Len-1)) * exp(-maxdist^2/(2*sigma_y^2));
%         J_y = (Len-psize_tmp)/sigma_y^2;
%         suminsep = 1 + sum(dist_y<(4*sigma_y)^2 & dist_y>(2*sigma_y)^2);
%         J_y = (Len-psize_tmp)/(sigma_y^2*suminsep);

        aff_y = exp(-dist_y./(2*sigma_y^2));
%         dist_y_from_x_F = aff_x_F'*aff_y;
        dist_y_from_x_F = corr_Vecs(aff_x_F,aff_y);
        
%         plot(aff_y(sind),aff(sind),'*'), axis equal, xlim([0 1]), xlim([0 1])
%         clf, plot(aff_x_F(sind),'*'), hold on, plot(aff(sind),'*')
        %%%%%%%%%%%%%%%%%%%%%%% Y to Y_PHI %%%%%%%%%%%%%%%%%%%%
%         st_sz = size(x_start,2); ph_sz = size(x_PHI,2);
%         yQ0_sz = size(y_from_Q0,2);
%         sz = min([st_sz,ph_sz,yQ0_sz])-1;
%         PHI_R =  x_start(:,1:sz)\x_PHI(:,1:sz);
%         PHI_T = x_PHI(:,end) - x_start(:,end);
%         
%         y_phi = zeros(F_2,sz+1);
% %         y_phi(:,end) = y_from_Q0(:,end) - PHI_T;
%         y_phi(:,1:end-1) = y_from_Q0(:,1:sz) / PHI_R;

        PHI_R = x_start\x_PHI;
        y_phi = y_from_Q0 / PHI_R;

        %%%%%%%%%%%%%%%%%%%%% Y_PHI to Y_F %%%%%%%%%%%%%%%%%%%%%%%

        theta = y_phi;
        iter_cnt = 0;
        while( iter_cnt < n_FLKOS_itrs )
            dist_y_F = feval(resfn,theta,X);
            [~,I] = sort(dist_y_F);
            inds = I( postBasedK_old-psize_tmp_old+1:postBasedK_old );
            psub = X(:,inds);
            iter_cnt = iter_cnt+1;
            if( iter_cnt < n_FLKOS_itrs )
                theta = feval(fitfn,psub);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Q(X_Start|Y_F) %%%%%%%%%%%%%%%%%%%%%%
        aff_y_F = exp(-dist_y_F./(2*sigma_old^2));
        aff_x_start = exp(-dist_old./(2*sigma_old^2));
%         dist_x_from_y_F = aff_y_F'*aff_x_start;
        dist_x_from_y_F = corr_Vecs(aff_y_F,aff_x_start);
        
        %%%%%%%%%%%%%%%%%%%%%%% MetropolisHastings %%%%%%%%%%%%%%%%%%%%
        MH_crit = (J_y/J_old)*dist_x_from_y_F/dist_y_from_x_F;
        MH_crit(isnan(MH_crit))=0;

        if  MH_crit > rand
            smpls_cnt = smpls_cnt + 1;
            Likli(smpls_cnt) = J_y;
            H_T(:,smpls_cnt) = aff_y;
            R_T(:,smpls_cnt) = dist_y;
            jmp_coin = mean(H_T(:,smpls_cnt));
            postBasedK_old = postBasedK;
            psize_tmp_old = psize_tmp;
            J_old = J_y;
            dist_old = dist_y;
            sigma_old = sigma_y;
            x_start = y_from_Q0;
        else
            Num_reject_4 = Num_reject_4 + 1;
        end
    else

        [~,I] = sort(dist_old);
        x_start_aff = exp(-dist_old/(2*sigma_old^2));
        pts_weights = x_start_aff;
        pts_weights(dist_old>(2*sigma_old)^2)=0;
        if max(pts_weights)==0
            pts_weights = x_start_aff;
        end
        
        Next_Mode_k = max([psize_tmp randsample(N, 1, true, eps+pts_weights(I))]);
        inds = I(Next_Mode_k - psize_tmp+1: Next_Mode_k);
        psub = X(:,inds);
%         for itr = 1:n_FLKOS_itrs
%             theta = feval(fitfn,psub);
%             dist = feval(resfn,theta,X);
%             [~,I] = sort(dist);
%             inds = I(k_m-psize+1:k_m);
%             psub = X(:,inds);
%         end
        x_PHI = feval(fitfn,psub);
        dist_x_PHI = feval(resfn,x_PHI,X);
        [~,I] = sort(dist_x_PHI);

        GI = G(I,I);
        Pts2HypAvgp = zeros(1,N);
        Pts2HypAvgp(1:k_m)=inf;
        subG1sum = sum(sum(GI(1:k_m,1:k_m)));
        subG2sum = sum(sum(GI(k_m+1:N,k_m+1:N)));
        cutG = sum(sum(GI(1:k_m,k_m+1:N))) + sum(sum(GI(k_m+1:N,1:k_m)));
        for kcnt = k_m + 1 : N-1
            vec = GI(:,kcnt);
            to_add_to_1 = 2*sum(vec(1:kcnt)) - vec(kcnt);
            to_take_from_2 = 2*sum(vec(kcnt:N)) - vec(kcnt);
            subG1sum = subG1sum + to_add_to_1;
            subG2sum = subG2sum - to_take_from_2;
            cutG = cutG - to_add_to_1 + to_take_from_2;
            cutCost = cutG/min([subG1sum subG2sum]);
            if cutCost>1.2*min(Pts2HypAvgp(1:kcnt-1))
                break;
            end
            Pts2HypAvgp(kcnt) = cutCost;
        end
        Pts2HypAvgp(1:k_m)=0;
        Pts2HypAvgp(isnan(Pts2HypAvgp))=0;
        Pts2HypAvgp(Pts2HypAvgp<0)=0;

        postBasedK = 12;
        if max(Pts2HypAvgp)>0
            postBasedK = randsample(N,1,true,Pts2HypAvgp);
        end
        postBasedK = min([postBasedK 40]);
        psize_tmp = max([psize floor(postBasedK/5)]);
        inds = I(postBasedK - psize_tmp+1: postBasedK);
        psub = X(:,inds);

        for itr = 1:n_FLKOS_itrs
            theta = feval(fitfn,psub);
            dist_y = feval(resfn,theta,X);
            [~,I] = sort(dist_y);
            inds = I(postBasedK-psize_tmp+1:postBasedK);
            psub = X(:,inds);
        end
        [sigma_y, ~, Len] = MSSE(dist_y, MSSEThresh, psize_tmp, postBasedK);
        if isinf(sigma_y)
            Num_reject_5 = Num_reject_5 + 1;    
            continue;
        end
%         maxdist = MaxSqRes^0.5;
        aff_y = exp(-dist_y./(2*sigma_y^2));
        J_y = exp(-mean(dist_y(I(1:Len)))/(2*(Len^0.5*sigma_y)^2))/(Len^0.5*sigma_y);
%         J_y = exp(-dist_y(I(postBasedK))/(2*sigma_y^2))/sigma_y;
%         J_y = (Len/(sigma_y*maxdist)) * ((normcdf(maxdist/sigma_y) - normcdf(-maxdist/sigma_y))^(Len-1)) * exp(-maxdist^2/(2*sigma_y^2));
%         J_y = (Len-psize_tmp)/sigma_y^2;
%         suminsep = 1 + sum(dist_y<(4*sigma_y)^2 & dist_y>(2*sigma_y)^2);
%         J_y = (Len-psize_tmp)/(sigma_y^2*suminsep);
        M_crit = (J_y/J_old);
%         clf, plot(dist_x_PHI(sind),'*'), hold on, plot(dist_y(sind),'*')
        
        smpls_cnt = smpls_cnt + 1;
        Likli(smpls_cnt) = J_y;
        H_T(:,smpls_cnt) = aff_y;
        R_T(:,smpls_cnt) = dist_y;
        jmp_coin = 0.2 + jmp_coin + mean(H_T(:,smpls_cnt));
        if  M_crit > rand
            postBasedK_old = postBasedK;
            psize_tmp_old = psize_tmp;
            J_old = J_y;
            dist_old = dist_y;
            sigma_old = sigma_y;
            x_start = theta;
        else
            Num_reject_7 = Num_reject_7 + 1;    
        end
    end
end
Likli = Likli(1:smpls_cnt);
H_T=H_T(:,1:smpls_cnt);
R_T=R_T(:,1:smpls_cnt);
% disp([Num_reject_1 Num_reject_2 Num_reject_3 Num_reject_4 Num_reject_5 Num_reject_6 Num_reject_7 smpls_cnt Total_jumps])

end

function corr = corr_Vecs(vec1, vec2)
    vec1 = vec1 - mean(vec1);
    vec1 = vec1 / std(vec1);
    vec2 = vec2 - mean(vec2);
    vec2 = vec2 / std(vec2);
    corr = mean(vec1.* vec2);
end