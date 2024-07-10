clear
close all
% clc

rng(0)

warning off

addpath('./model_specific');
addpath('./algorithmSubFunctions');

fundlist = dir('./data/adelaidermf');
fundLabel = {};
for cnt = 3 : length(fundlist)
    fundLabel(cnt-2) = mat2cell([fundlist(3).folder '/' fundlist(3).name], 1);
end

numRun = 100;
results = zeros(length(fundLabel), 4);
for seq_num = 1:length(fundLabel)
    seq_fpath = cell2mat(fundLabel(seq_num));
    disp(['running seq: ', seq_fpath])
    load(seq_fpath);
    numModels = max(label) - min(label) + 1;

    numPoints=zeros(max(label), 1);
    cnt = 1;
    for i=min(label):max(label)
       numPoints(cnt) = sum(label==i);
       cnt = cnt + 1;
    end
    disp('Number of points(outliers first): ')
    disp(numPoints)

    %Parameter Declaration
    model_type = 'fundamental';

    %remove repeating rows in data
    [data,ia,ic] = unique(data','rows');
    data = data';
    label = label(ia);

    dat_img_1 = normalise2dpts(data(1:3,:));
    dat_img_2 = normalise2dpts(data(4:6,:));

    X = [dat_img_1; dat_img_2];

    n_once_double_triple_list = [2, 4, 8];  %choose 2 for once, 4 for double and 8 for triple
    n_once_double_triple = n_once_double_triple_list(1);

    numHypo = floor(length(unique(label)) * n_once_double_triple * log(length(X)));
    %50 times less than what is written in the paper because that is for
    %line fitting which is faster and maybe necessary.

    miss_rateH = zeros(1,numRun);
    ttimeH = zeros(1,numRun);
    for nRun=1:numRun
        
        tic;
        if(0) %IPTA'16
            [H,R] = MCMC_km_funda(X, model_type, label, numHypo);
            ClustLabels = SC_with_outliers_funda(H, R, numModels, model_type, X, label);
        else  %CVIU'17
            NUMBAYESITR = 5;
            G = nan;
            Liks_MAT = [];  
            H_MAT = [];
            R_MAT = [];
            for Bayesitr = 1:NUMBAYESITR
                
                %MCNC uses FLKOS and samples k from NCut cost liklihood
                [Hp,Rp,Liks] = MCNC_funda(X, model_type, G, numHypo, label);

                %SLKOS is my idea to allow FLKOS handle skewed densities
                % I propose MSNC: M_CMC S_LKOS N_ormalized C_ut
%                 [Hp,Rp,Liks] = MSNC_funda(X, model_type, G, numHypo, label);
                
%                 [Hp,Rp,Liks] = PerfectSamples(Tracks, model_type, G, n_samples,s);
    
                Liks_MAT = [Liks_MAT Liks];
                H_MAT = [H_MAT Hp];
                R_MAT = [R_MAT Rp];
    
                Dv = sum(H_MAT*diag(Liks_MAT),2);
                De = sum(H_MAT,1);
                DvIh = diag(Dv.^-0.5);
                DeI = diag(De.^-1);
                G = DvIh*H_MAT*DeI*H_MAT'*DvIh;
    
                Dv_OUT = Dv<max(Dv)*0.02;
                G(Dv_OUT,:)=0;
                G(:,Dv_OUT)=0;
                G = G/max(G(:));            
                G(isnan(G))=0;
            end
            ClustLabels = SC_including_outliers_28( G, numModels, model_type, X, label');
        end
        ttime = toc;

        miss_rate = 100 - compute_accuracy(label+1, ClustLabels');

        miss_rateH(nRun) = miss_rate;
        ttimeH(nRun) = ttime;
        disp(['misclass error = ', num2str(miss_rate)])
        disp(['Time = ', num2str(ttimeH(nRun))])
    end
    myresults(:,seq_num) = miss_rateH;
    results(seq_num,:) = [
        mean(miss_rateH), median(miss_rateH), mean(ttimeH), median(ttimeH)];
end

%%
function vis_sample(img1, data, label, ClustLabels, miss_rate)
    figure
    subplot 121
    imshow(img1);hold on
    gscatter(data(1,:), data(2,:), label,[],[],20)
    title('True Clusters')
    subplot 122
    imshow(img1);hold on
    gscatter(data(1,:), data(2,:), ClustLabels,[],[],20)
    title(['Estimated Clusters, acc: ' num2str(miss_rate)])
end