function [ ClustLabels, Energy] = SC_including_outliers_28( G, numModels, model_type, X, label)
% [~,sind]=sort(label); % use when visualizing by labels
[fitfn, resfn,~,psize] = getModelParam(model_type);
psize = psize + 2;
N = size(G,1);

% 1 - Each eigenvector includes two structures, one on the positive side
% and one on the negative. We keep the one with larger components, if
% both have same maximum we keep them both as seperate structures.
% 2 - Bridges cannot apear later, no eigenvec can bridge
% between eigvecs with larger eignevalues
% 3 - We rather have oversegmentation than under. so we get rid
% of stronger bridges as well so if the first eigenvec is a bridge
% for graph connectivity it will go away.
% this leaves us with oversegments that we megre to reach the 
% desired number of strucutres
% 4 - Those data points that are not expressed are outliers. if
% we have prior info that there are no outliers, we add them one by one
% to closest structures and remodel the structures until no data point is
% left out.
% 5 - Running KMeans is absoulutly meaningless. We will be just looking at
% how close data points are to eigenvectors. The truth is that the 
% components are widely spread along each axes. But there mught be huge
% gaps between them on an axis or many axes or points may be 
% concentraited close to origin. Kmeans won't work in this situation.
% The solution is to see which eigenvector, the data point belongs to by: on 
% which it has its maximum coefficient. End of story. End of my PhD.
% I wasted so much time on the concept of Kmeans.


THR_HighProb = normpdf(1,0,1)/normpdf(0,0,1);     % from 1 to 2 sigma is for inliers
THR_inlier = normpdf(2,0,1)/normpdf(0,0,1);     % from 1 to 2 sigma is for inliers
THR_GrossOut = normpdf(3,0,1)/normpdf(0,0,1);   % from 2.75 is for gross outliers


[U,sigs,~]=svd(G);
sigs = diag(sigs);
sigs = sigs/sigs(1);
sigs_cumsum = cumsum(sigs);
d = find(sigs_cumsum>numModels,1);
d = d - 1;

%     Uz = U(:,1:numModels);
%     [~, ClustLabels] = vl_kmeans(Uz', numModels, 'distance', 'l1', ...
%      'algorithm', 'elkan', 'NumRepetitions', 200, 'NumTrees', 3); 
%     ClustLabels = double(ClustLabels');
% return

Uz = zeros(N,1);
maxG = max(G,[],2);
maxG = maxG/max(maxG);
vec2Add = 0;
done = 0;
if d==numModels
    segcnt = 0;    
else
    segcnt = 1;
end
while ~done
    segcnt = segcnt + 1;
    eigVecPos = U(:,segcnt);
    eigVecPos = eigVecPos/max(eigVecPos);
    eigVecPos(eigVecPos<0)=0;
    eigVecPos(isnan(eigVecPos))=0;
    
    eigVecNeg = -U(:,segcnt);
    eigVecNeg = eigVecNeg/max(eigVecNeg);
    eigVecNeg(eigVecNeg<0)=0;
    eigVecNeg(isnan(eigVecNeg))=0;
    
    if max(eigVecPos)>0.5*max(eigVecNeg)
        if sum(eigVecPos>THR_inlier)>psize
            vec2Add = vec2Add + 1;
            Uz(:,vec2Add)=eigVecPos/max(eigVecPos);
        end
    end
    
    if max(eigVecNeg)>0.5*max(eigVecPos)
        if sum(eigVecNeg>THR_inlier)>psize
            vec2Add = vec2Add + 1;
            Uz(:,vec2Add)=eigVecNeg/max(eigVecNeg);
        end
    end
    
    maxUz = max(Uz,[],2);
    maxUz = maxUz/max(maxUz); maxUz(isnan(maxUz))=0;
    try
        if ( ( (segcnt>=numModels) && (min(maxUz(maxG>0.01))>0.1) ) || ...
             (segcnt>=max([numModels min([2*numModels-1  d])])) )
            done = 1;
        end
    catch
        disp('here')
    end
end
Uz(Uz<THR_inlier)=0;

d = size(Uz,2);
if d>numModels
    Mrec = zeros(1,d);
    for segcnt = 3 : d
        Uzp = Uz(:,1:segcnt);
        Uzp(Uzp>0.1)=1;
        Uzp(Uzp<=0.1)=0;
        OvLp = Uzp'*Uzp;
        OvLp = OvLp./repmat(sum(Uzp),segcnt,1);
        OvLp(OvLp<0.01)=0;

        M1 = OvLp(1:end-1,1:end-1);
        M2 = (OvLp(end,1:end-1)'*OvLp(end,1:end-1)).^0.5;
        M = M2./M1;
        M(isnan(M))=0;
        if ~isempty(M)
            Mrec(segcnt) = max(M(:));
        end
    end
    if sum(Mrec<4)<numModels
        id = find(Mrec>=4);
        Mrec(id(1:numModels-sum(Mrec<4)))=0;
    end
    Uz = Uz(:,Mrec<4);
    d = size(Uz,2);
end
if d>numModels
    Mrec = zeros(1,d);
    for segcnt = 1 : d-2
        Uzp = Uz(:,segcnt:end);
        Uzp(Uzp>0.1)=1;
        Uzp(Uzp<=0.1)=0;
        OvLp = Uzp'*Uzp;
        OvLp = OvLp./repmat(sum(Uzp),d-segcnt+1,1);
        OvLp(OvLp<0.01)=0;

        M1 = OvLp(2:end,2:end);
        M2 = (OvLp(1,2:end)'*OvLp(1,2:end)).^0.5;
        M = M2./M1;
        M(isnan(M))=0;
        if ~isempty(M)
            Mrec(segcnt) = max(M(:));
        end
    end
    
    if sum(Mrec<4)<numModels
        id = find(Mrec>=4);
        Mrec(id(1:numModels-sum(Mrec<4)))=0;
    end

    Uz = Uz(:,Mrec<4);
    d = size(Uz,2);
end

Uz(:,end+1)=THR_GrossOut;
[~,ClustLabels] = max(Uz,[],2);
Uz = Uz(:,1:end-1);

Uzp = Uz;
for clcnt = 1 : d
    inds = ClustLabels==clcnt;
    [~,I] = sort(Uzp(:,clcnt),'descend');
    if sum(inds)<psize
        ClustLabels(I(1:psize))=clcnt;
    end
    Uzp(I(1:psize),:)=0;
end

%%%%%%%%%%%% Rank Reduction %%%%%%%%%%%
if d>numModels
    sqScaledDist = zeros(N,d);
    for segcnt = 1 : d
        clus_inds = ClustLabels==segcnt;
        psub = X(:,clus_inds);
        theta = feval(fitfn,psub);
        sqDist = feval(resfn,theta,X);
        outvar = sum(sqDist(clus_inds))/(sum(clus_inds)-psize);
        sqScaledDist(:,segcnt) = (sqDist/outvar).^0.5;
    end
end

while(d>numModels)
    
    Weighted_Dist = inf + zeros(d);
    Gbased = zeros(d);
    for segcnt_i = 1 : d
        for segcnt_j = 1 : d
            Weighted_Dist(segcnt_i, segcnt_j) = mean(sqScaledDist(ClustLabels==segcnt_j, segcnt_i));
            Gbased(segcnt_i, segcnt_j) =  mean(mean(G(ClustLabels==segcnt_i,ClustLabels==segcnt_j)));
        end
    end
    Gbased = Gbased - diag(diag(Gbased));
    Gw = 2./(Weighted_Dist + Weighted_Dist');
    if max(Gbased(:))>0
        Gw = Gbased.^0.5.*Gw;
    else
        Gw = Gw - diag(diag(Gw));
    end

    
    maxGz = max(Gw(:));
    [idi,idj] = find(Gw==maxGz);
    idi = idi(1); idj = idj(1);

    idlow = min([idi idj]);
    idhigh = max([idi idj]);
%     subplot(1,2,1), plot(ClustLabels(sind),'*'), subplot(1,2,2), imagesc(Gw)
    ClustLabels(ClustLabels==idhigh) = idlow;
    for clcnt = idhigh + 1 : max(ClustLabels)
        ClustLabels(ClustLabels==clcnt)=clcnt-1;
    end
    
    clus_inds = ClustLabels==idlow;
    psub = X(:,clus_inds);
    theta = feval(fitfn,psub);
    sqDist = feval(resfn,theta,X);
    outvar = sum(sqDist(clus_inds))/(sum(clus_inds)-psize);
    sqScaledDist(:,idlow) = (sqDist/outvar).^0.5;

    sqScaledDist(:,idhigh)=[];
    d = d - 1;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indx = ClustLabels == d + 1;
if sum(indx)>0
    newU = zeros(N,d);
    newR = zeros(N,d);
    for cnt = 1 : d
        clus_inds = ClustLabels==cnt;
        psub = X(:,clus_inds);
        theta = feval(fitfn,psub);
        sqDist = feval(resfn,theta,X);
        outvar = 2*median(sqDist(clus_inds));
        newU(:,cnt) = sqDist/outvar;
        newR(:,cnt) = sqDist;
    end
    
    while sum(indx)>0
        Clslbl = ClustLabels;
        minSqScDist = inf + zeros(N,1);
        [minSqScDist(indx), Clslbl(indx)] = min(newU(indx,:),[],2);
        [~,ptid] = min(minSqScDist);
        ptid = ptid(1);
        clsnum = Clslbl(ptid);
        ClustLabels(ptid) = clsnum;
%         plot(ClustLabels(sind),'*'), pause(.1)
        inds = ClustLabels == clsnum;
        psub = X(:,inds);
        theta = feval(fitfn,psub);
        sqDist = feval(resfn,theta,X);
        outvar = 2*median(sqDist(inds));
        newU(:,clsnum) = sqDist/outvar;
        newR(:,clsnum) = sqDist;
        indx = ClustLabels == d + 1;
    end
    Energy = mean(newR(:)); 
    % figure, plot(ClustLabels(sind),'*')

end

end