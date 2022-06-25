%% Data treatment
data = readtable('agaricus-lepiota.txt','ReadVariableNames',false);
cats = categorical(data{:,:});

%remove column with '?'
data_clean= removevars(data,{'Var12'});

% check if some '?' is still there
ismember('?',data_clean{:,:})

% categorize features and classes
data_cats = categorical(data_clean{:,:});
labels_cat = data_cats(:,1);
features_cat = data_cats(:,2:width(data_cats));

%% Categorical to double
labels = grp2idx(labels_cat);
features = zeros(size(features_cat));
for i = 1:width(features_cat)
    features(:,i) = grp2idx(features_cat(:,i));
end

%% Similarity Matrix
similarity_matrix_hamming = pdist2(features, features, 'cosine');

%% Threshold
hamming_thr = 0.1;

similarity_matrix_hamming(similarity_matrix_hamming < hamming_thr) = 0.1;

%% D matrix, D^(-1) and D^(-1/2)
D_hamming = zeros(height(features));
D_hamming_invsqr = zeros(height(features));
D_hamming_inv = zeros(height(features));
for i = 1:height(features)
    D_hamming(i,i) = sum(similarity_matrix_hamming(i,:)>0);
    D_hamming_invsqr(i,i) = 1/sqrt(D_hamming(i,i));
    D_hamming_inv(i,i) = 1/D_hamming(i,i);
end

%% Eigenvalues and eigenvectores for the type of Laplacian used
%1 for Laplacian symmetric normalized
%2 for Laplacian random walk
%3 for Laplacian unnormalized
laplacian_type = 2;

switch laplacian_type
    case 1
        Lsn_hamming = D_hamming_invsqr * L_hamming * D_hamming_invsqr;
        [EigVec_hamming, EigVal_hamming] = eigs(Lsn_hamming,2,'smallestabs');
    case 2
        Lrw_hamming = D_hamming_inv * L_hamming;
        [EigVec_hamming, EigVal_hamming] = eigs(Lrw_hamming,2,'smallestabs');
    case 3
        L_hamming = D_hamming - similarity_matrix_hamming;
        [EigVec_hamming, EigVal_hamming] = eigs(L_hamming,2,'smallestabs');
end
%% Using k-means
k=2;
idx_hamming = kmeans(real(EigVec_hamming(:,1:2)), k);

%%
idx_hamming(idx_hamming == 1) = 14;
idx_hamming(idx_hamming == 2) = 5;

%%
ok = idx_hamming == labels;
sum(ok)

%%
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

%% Performance
TP = CM_knn_H(1,1);
FP = CM_knn_H(1,2);
FN = CM_knn_H(2,1);
TN = CM_knn_H(2,2);
perf.error = (FP+FN)/(TP+FP+TN+FN);
perf.accuracy = (TP+TN)/(TP+FP+TN+FN);
perf.precision = TP/(TP+FP);
perf.recall = TP/(TP+FN);
perf.specificity = TN/(TN+FP);
perf.f_score = 2*((precision*recall)/(precision+recall));

P_isClass1_ClassifiedAsClass1 = TP/(TP+FP)
P_isClass2_ClassifiedAsCLass1 = FP/(TP+FP)
P_isClass1_ClassifiedAsClass2 = FN/(FN+TN)
P_isClass2_ClassifiedAsClass2 = TN/(FN+TN)