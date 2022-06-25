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

%% Similarity Matrix and threshold
metric = 'hamming'

similarity_matrix = pdist2(features, features, metric);

hamming_thr = 0.4;
cosine_thr = 0.15;

if strcmp(metric, 'hamming')
    similarity_matrix = 1 - similarity_matrix;
    similarity_matrix = similarity_matrix - min(min(similarity_matrix));
    similarity_matrix(similarity_matrix < hamming_thr) = 0.0;
end

if strcmp(metric, 'cosine')
    similarity_matrix(similarity_matrix < cosine_thr) = 0.0;
end

%% D matrix, D^(-1) and D^(-1/2)
D = zeros(height(features));
D_invsqr = zeros(height(features));
D_inv = zeros(height(features));
for i = 1:height(features)
    D(i,i) = sum(similarity_matrix(i,:)>0);
    D_invsqr(i,i) = 1/sqrt(D(i,i));
    D_inv(i,i) = 1/D(i,i);
end

%% Eigenvalues and eigenvectores for the type of Laplacian used
%1 for Laplacian symmetric normalized
%2 for Laplacian random walk
%3 for Laplacian unnormalized
laplacian_type = 1;

L = D - similarity_matrix;

switch laplacian_type
    case 1
        Lsn = D_invsqr * L * D_invsqr;
        [EigVec, EigVal] = eigs(Lsn,2,'smallestabs');
    case 2
        Lrw = D_inv * L;
        [EigVec, EigVal] = eigs(Lrw,2,'smallestabs');
    case 3
        [EigVec, EigVal] = eigs(L,2,'smallestabs');
end
%% Using k-means
k=2;
idx = kmeans(real(EigVec(:,1:2)), k);

%%
idx(idx == 1) = 14;
idx(idx == 2) = 5;

%%
ok = idx == labels;
sum(ok)

%%
figure
CM_knn_H=confusionmat(labels,idx)
confusionchart(CM_knn_H,'FontSize',20)

%%
figure
colormap winter

subplot(2,1,1)
scatter(EigVec(:,1),EigVec(:,2), 4, labels, 'filled')
subplot(2,1,2)
scatter(EigVec(:,1),EigVec(:,2), 4, idx, 'filled')

%% Performance
TP = CM_knn_H(1,1);
FP = CM_knn_H(1,2);
FN = CM_knn_H(2,1);
TN = CM_knn_H(2,2);

P_isClass1_ClassifiedAsClass1 = TP/(TP+FP)
P_isClass2_ClassifiedAsCLass1 = FP/(TP+FP)
P_isClass1_ClassifiedAsClass2 = FN/(FN+TN)
P_isClass2_ClassifiedAsClass2 = TN/(FN+TN)

error = (FP+FN)/(TP+FP+TN+FN)
accuracy = (TP+TN)/(TP+FP+TN+FN)
precision = TP/(TP+FP)
recall = TP/(TP+FN)
specificity = TN/(TN+FP)
f_score = 2*((precision*recall)/(precision+recall))