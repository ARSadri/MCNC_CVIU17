function [ ClustLabels, method_id ] = Cut_Aff( H, numModels, DBSCAN_Effort )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
numModels = numModels + 1;

H(H<normpdf(5,0,1))=0;
G = corr_H2G(H); G(G<0)=0;
[V,~,~]=svd(G);
V = V(:,1: numModels);

ClustLabels = 0;
ClLabels = zeros(size(H,1),DBSCAN_Effort);
th = zeros(DBSCAN_Effort,1);
while max(ClustLabels) ~= numModels && DBSCAN_Effort>0
    th(DBSCAN_Effort) = max(V(:))*rand/2;
    [~,ClustLabels] = dbscan(V', th(DBSCAN_Effort), 6);
    ClLabels(:,DBSCAN_Effort)=ClustLabels;
    DBSCAN_Effort = DBSCAN_Effort - 1;
end

m_ClLabels = max(ClLabels);
m_ClLabels(m_ClLabels<numModels)=inf;
[numclsdb,id]= min(m_ClLabels);

method_id = 1;
if isinf(numclsdb)
    disp('DBSCAN failed')
    [~, ClustLabels, ~] = vl_kmeans(V', numModels, 'distance', 'l1', ...
         'algorithm', 'elkan', 'NumRepetitions', 100, 'NumTrees', 3); 
    ClustLabels = double(ClustLabels');
else
    ClustLabels = ClLabels(:,id);
    
    if numclsdb > numModels
    
    % merging
        numclss = max(ClustLabels);
        moment = zeros(numclss,numModels);
        for segcnti = 1 : numclss - 1
            moment(segcnti,:) = mean(V(ClustLabels==segcnti,:));
        end

        [~, segLabels, ~] = vl_kmeans(moment', numModels, 'distance', 'l1', ...
             'algorithm', 'elkan', 'NumRepetitions', 100, 'NumTrees', 3); 
        segLabels = double(segLabels');

        for segcnter = 1 : length(segLabels)
            ClustLabels(ClustLabels==segcnter)=-segLabels(segcnter);
        end
        ClustLabels = -ClustLabels;
        method_id = 0;
    end
end

end

