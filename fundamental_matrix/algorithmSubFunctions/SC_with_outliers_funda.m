function [ outLabels ] = SC_with_outliers_funda( H, R, numModels, model_type, X, label, G)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

INSPTHR = 2;
SVD_CUT_THR = 0.03;
[ fitfn, resfn,~,psize] = getModelParam(model_type);

N = size(H,1);
HS = sort(H','descend');
H(H<normpdf(2,0,1))=0;
min_num_hyps = ceil(size(H,2)/99);
crit = mean(HS(1:min_num_hyps,:));   %1% of hyps
inds = crit >= mean(crit(crit<median(crit)));
if isnan(G)
    Np = sum(inds);
    H = H(inds,:);
    R = R(inds,:);
    GH = corr_H2G(H); GH(GH<0)=0;
    RH = corr_H2G(R); RH(RH<0)=0; RH(RH>0)=1;
    G = GH.* RH;
end

[U,sV,~]=svd(G);
sV = diag(sV)';
ss = sV./[1 cumsum(sV(1:end-1))];
[~, d] = max( ss'<SVD_CUT_THR);
d = max(d,numModels);
% d = numModels * 2;
Ur = U(:,1: d);
[~, ClustLabels] = vl_kmeans(Ur', d, 'distance', 'l1', ...
     'algorithm', 'elkan', 'NumRepetitions', 200, 'NumTrees', 3); 
ClustLabels = double(ClustLabels');

%% visualization of over-segmentation
% a = (1 : length(label))';
% b=[];
% for segcnt = min(label):max(label)
%     b = [b;a(label==segcnt)];
% end
% mylabel = label(inds==1);
% mya = (1 : length(mylabel))';
% myb=[];
% for segcnt = min(mylabel):max(mylabel)
%     myb = [myb;mya(mylabel==segcnt)];
% end
% cllbls_before = zeros(N,1);
% cllbls_before(inds) = ClustLabels;
% figure, plot(cllbls_before(b),'*'), xlim([0 N])

ClustLabels_copy = ClustLabels;

Xp = X(:,inds);
ClustLabels = ClustLabels(inds);
d = max(ClustLabels);
while(d>numModels)
    spr_clsn = zeros(d,d);
    dist_clsn = zeros(d,d);
    for segcnti = 1 : d
        indx = ClustLabels==segcnti;
        psub = Xp(:,indx);
        theta = feval(fitfn,psub);
        Sdist = feval(resfn,theta,Xp);
        sigma = (sum(Sdist(indx))/(sum(indx)-psize))^0.5;
        for segcntj = 1 : d
            Sdisttmp = Sdist(ClustLabels==segcntj).^0.5;
            Sdisttmp_Sorted = sort(Sdisttmp);
            dist_cls = Sdisttmp_Sorted(ceil(0.5*end));
            dist_clsn(segcnti,segcntj) = dist_cls;
            spr_clsn(segcnti,segcntj) = dist_cls/sigma;
        end
    end
    clear d_clsp
    d_clsp(:,:,1) = spr_clsn; d_clsp(:,:,2) = spr_clsn'; d_clsp = max(d_clsp,[],3);
    d_clsp = d_clsp + 10^10*eye(d);
    inspcrit = min(d_clsp(:));
    if inspcrit<=INSPTHR
        inspcrit
        [idi,idj]=find(d_clsp==inspcrit,1);
    else
        clear d_clsp
        d_clsp(:,:,1) = dist_clsn; d_clsp(:,:,2) = dist_clsn'; d_clsp = max(d_clsp,[],3);
        d_clsp = d_clsp + 10^10*eye(d);
        
        dist_clsnz = dist_clsn + 10^10*eye(d);
        dist_clsnp = dist_clsn;
        dist_clsnp = dist_clsnp - diag(diag(dist_clsnp));
        Hypsenrg = zeros(d);
        for clcnti = 1 : d
            for clcntj = 1 : d
                the_row = min([dist_clsnz(clcnti,:);dist_clsnz(clcntj,:)]);
                sum_r = sum(the_row) - the_row(clcnti) - the_row(clcntj);
                the_clm = min( [dist_clsnz(:,clcnti), dist_clsnz(:,clcntj)] ,[],2);
                sum_cl = sum(the_clm) - the_clm(clcnti) - the_clm(clcntj);
                dists = dist_clsnp;
                dists([clcnti;clcntj],:)=[];
                dists(:,[clcnti;clcntj])=[];
                mysum = sum_r + sum_cl + sum(dists(:));
                Hypsenrg(clcnti,clcntj)=mysum/(d^2-3*d);
            end
        end
        d_clsz = d_clsp - Hypsenrg;
        [idi,idj]=find(d_clsz==min(d_clsz(:)),1);
    end

    ClustLabels(ClustLabels==idj) = idi;
    ClustLabelsp = 0*ClustLabels;
    cnt = 0;
    for clcnt = unique(ClustLabels)'
        cnt = cnt + 1;
        ClustLabelsp(ClustLabels==clcnt)=cnt;
    end
    ClustLabels = ClustLabelsp;

    d = max(ClustLabels);
    
%     cllbls_before = zeros(N,1);
%     cllbls_before(inds) = ClustLabels;
%     figure, plot(cllbls_before(b),'.'), xlim([0 N])
%     getframe
end

cllbls = zeros(N,1);
cllbls(inds) = ClustLabels;
if max(cllbls) == 0
    disp('asdf')
end
Hout = zeros(N, max(cllbls));
for segcnt = 1 : max(cllbls)
    indx = cllbls==segcnt;
    psub = X(:,indx);
    theta = feval(fitfn,psub);
    Sdist = feval(resfn,theta,X);
    sigma = (sum(Sdist(indx))/(sum(indx)-psize))^0.5;
    Hout(:,segcnt) = exp(-Sdist/(2*sigma^2));
end
[~,outLabels] = max(Hout,[],2);
cllbls(cllbls==0)=outLabels(cllbls==0);
outLabels = cllbls;

end