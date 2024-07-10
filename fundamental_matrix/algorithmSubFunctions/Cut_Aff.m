function [ ClustLabels ] = Cut_Aff( H, numModels, label )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
numModels = numModels + 1;

H(H<normpdf(5,0,1))=0;
G = corr_H2G(H); G(G<0)=0;
[V,~,~]=svd(G);
V = V(:,1: numModels);


ClustLabels = 0;
while max(ClustLabels) < numModels
    th = max(V(:))*rand/10;
    [~,ClustLabels] = dbscan(V', th, 6);
end

numclss = max(ClustLabels);

if numclss > numModels     %merging
    
    moment = zeros(numclss,numModels);
    for segcnti = 1 : numclss
        moment(segcnti,:) = mean(V(ClustLabels==segcnti,:));
    end

%     cldist = zeros(numclss,numclss);
%     for segcnti = 1 : numclss - 1
%         for segcntj = segcnti + 1 : numclss
%             data_i = V(ClustLabels==segcnti,:);
%             data_j = V(ClustLabels==segcntj,:);
%             dists_i_to_j = sum((data_i - repmat(moment(segcntj,:),size(data_i,1),1)).^2,2);
%             [~,closest_pt_from_i_to_j] = min(dists_i_to_j);
%             dists_j_to_i = sum((data_j - repmat(moment(segcnti,:),size(data_j,1),1)).^2,2);
%             [~,closest_pt_from_j_to_i] = min(dists_j_to_i);
%             cldist(segcnti,segcntj) = sum((data_i(closest_pt_from_i_to_j,:) - data_j(closest_pt_from_j_to_i,:)).^2)^0.5;
%         end
%     end
%     cldist = cldist+ cldist';
    
    [~, segLabels, ~] = vl_kmeans(moment', numModels, 'distance', 'l1', ...
         'algorithm', 'elkan', 'NumRepetitions', 100, 'NumTrees', 3); 
    segLabels = double(segLabels');

    for segcnter = 1 : length(segLabels)
        ClustLabels(ClustLabels==segcnter)=-segLabels(segcnter);
    end
    ClustLabels = -ClustLabels;
end