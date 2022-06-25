data = readtable('agaricus-lepiota.txt','ReadVariableNames',false);
cats = categorical(data{:,:});

%remove column with '?'
data_clean= removevars(data,{'Var12'});

% check if some '?' is still there
ismember('?',data_clean{:,:});

% categorize features and 
data_cats = categorical(data_clean{:,:});
labels_cat = data_cats(:,1);
features_cat = data_cats(:,2:width(data_cats));
%%
labels = grp2idx(labels_cat);
features = zeros(size(features_cat));
for i = 1:width(features_cat)
    features(:,i) = grp2idx(features_cat(:,i));
end

%%
similarity_matrix_hamming = pdist2(features, features, 'hamming');
%%
hamming_thr = 0.5;

similarity_matrix_hamming(similarity_matrix_hamming < hamming_thr) = 0.0;

%%
D_hamming = zeros(height(features));
D_hamming_invsqr = zeros(height(features));
D_hamming_inv = zeros(height(features));
for i = 1:height(features)
    D_hamming(i,i) = sum(similarity_matrix_hamming(i,:)>0);
    D_hamming_invsqr(i,i) = 1/sqrt(D_hamming(i,i));
    D_hamming_inv(i,i) = 1/D_hamming(i,i);
end
%%
L_hamming = D_hamming - similarity_matrix_hamming;

Lsn_hamming = D_hamming_invsqr * L_hamming * D_hamming_invsqr;

Lrw_hamming = D_hamming_inv * L_hamming;

%% Laplacian symmetric normalized
k = 2;
[EigVec_hamming, EigVal_hamming] = eig(Lsn_hamming);

%% Laplacian random walk
k = 2;
[EigVec_hamming, EigVal_hamming] = eig(Lrw_hamming);
%% Laplacian unnormalized
k = 2;
[EigVec_euclidean, EigVal_euclidean] = eig(L_euclidean);
[EigVec_hamming, EigVal_hamming] = eig(L_hamming);

%% sort eigenvalues and get 2 smaller
[D,I] = sort(diag(EigVal_hamming));
EigVec_hamming = EigVec_hamming(:, I);

[D,I] = sort(diag(EigVal_euclidean));
EigVec_euclidean = EigVec_euclidean(:, I);
%%
idx_euclidean = kmeans(EigVec_euclidean(:,1:2), k);
idx_hamming = kmeans(real(EigVec_hamming(:,1:2)), k);
%%
idx_euclidean(idx_euclidean == 1) = 14;
idx_euclidean(idx_euclidean == 2) = 5;

idx_hamming(idx_hamming == 1) = 14;
idx_hamming(idx_hamming == 2) = 5;

%%
ok = idx_hamming == labels;
sum(ok)

ok = idx_euclidean == labels;
sum(ok)

%%
figure
CM_knn_E=confusionmat(labels,idx_euclidean)
confusionchart(CM_knn_E,'FontSize',20)

figure
CM_knn_H=confusionmat(labels,idx_hamming)
confusionchart(CM_knn_H,'FontSize',20)

%%
figure
colormap winter

subplot(2,1,1)
scatter(EigVec_hamming(:,1),EigVec_hamming(:,2), 4, labels, 'filled')
subplot(2,1,2)
scatter(EigVec_hamming(:,1),EigVec_hamming(:,2), 4, idx_hamming, 'filled')

figure
colormap winter

subplot(2,1,1)
scatter(EigVec_euclidean(:,1),EigVec_euclidean(:,2), 4, labels, 'filled')
subplot(2,1,2)
scatter(EigVec_euclidean(:,1),EigVec_euclidean(:,2), 4, idx_euclidean, 'filled')

%%