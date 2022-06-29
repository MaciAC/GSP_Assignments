%% 1. Reading and cleaning up the data
close all
clear all
clc

T = readtable('LazegaLawyers/ELwork36.dat');

%% 2. Plot the graph
A = table2array(T);
G = graph(A);
figure(1);
plot(G,'layout','force');

%% 3. Common Neighbors

N_nodes = height(T);
N_edges = N_nodes * (N_nodes - 1) / 2;
Adjacency_orig = zeros(N_nodes);

%create list of edges
edges = zeros(N_edges, 2);
selected_idx = randperm(N_edges, N_edges/5); 
n = 1;
for i = 1:N_nodes
    for j = i+1:N_nodes
        edges(n,:) = [i j];
        n=n+1;
    end
end

score_function = 'Common_Neighbors';
score_function = 'Jaccard_index';
score_function = 'Addamic_Addar';
for i = selected_idx    
    switch score_function
        case 'Common_Neighbors'
            thr_ = 1:10;
            Adjacency_orig(edges(i,1), edges(i,2)) = Common_Neighbors(A(edges(i,1),:), A(edges(i,2),:));
        case 'Jaccard_index'
            thr_ = 0:0.1:1;
            Adjacency_orig(edges(i,1), edges(i,2)) = Jaccard_index(A(edges(i,1),:), A(edges(i,2),:));
        case 'Addamic_Addar'
            thr_ = 0:0.2:2;
            Adjacency_orig(edges(i,1), edges(i,2)) = Addamic_Addar(A(edges(i,1),:), A(edges(i,2),:));
            Adjacency_orig(Adjacency_orig == Inf) = -Inf;
            Adjacency_orig(Adjacency_orig == -Inf) = max(max(Adjacency_orig));
    end
    Adjacency_orig(edges(i,2), edges(i,1)) = Adjacency_orig(edges(i,1), edges(i,2));
end

%% 4. Probability of detection and false alarm

A_roc = zeros(width(selected_idx),10);
Adjacency_roc = zeros(width(selected_idx),10);
P_detection = zeros(10,1);
P_false_alarm = zeros(10,1);

k=1;
for thr = thr_

Adjacency = Adjacency_orig;
Adjacency(Adjacency_orig < thr) = 0.0;


Detection = Adjacency .* A;
temp = zeros(N_nodes);
temp(Adjacency > 0) = 1;
False_Alarm = A - temp;
n_detected = 0;
n_expected = 0;
zero_expected = 0;
n_false_alarm = 0;
j=1;

for i = selected_idx
    A_roc(j,k) =  A(edges(i,1), edges(i,2));
    Adjacency_roc(j,k) = Adjacency(edges(i,1), edges(i,2));
    if A(edges(i,1), edges(i,2)) > 0
        n_expected = n_expected + 1;
    else
        zero_expected = zero_expected + 1;
    end
    if Detection(edges(i,1), edges(i,2)) > 0
        n_detected = n_detected + 1; 
    end
    if False_Alarm(edges(i,1), edges(i,2)) == -1
        n_false_alarm = n_false_alarm + 1;
    end
    j = j +1;
end

P_detection(k,1) = n_detected / n_expected;

P_false_alarm(k,1) = n_false_alarm / zero_expected;
k = k + 1;
end

figure
plotroc(A_roc(:,1)', Adjacency_roc(:,1)','Threshold 1',A_roc(:,2)', Adjacency_roc(:,2)','Threshold 2',...
    A_roc(:,3)', Adjacency_roc(:,3)','Threshold 3',A_roc(:,4)', Adjacency_roc(:,4)','Threshold 4',...
    A_roc(:,5)', Adjacency_roc(:,5)','Threshold 5',A_roc(:,6)', Adjacency_roc(:,6)','Threshold 6',...
    A_roc(:,7)', Adjacency_roc(:,7)','Threshold 7',A_roc(:,8)', Adjacency_roc(:,8)','Threshold 8',...
    A_roc(:,9)', Adjacency_roc(:,9)','Threshold 9',A_roc(:,10)', Adjacency_roc(:,10)','Threshold 10');

%% functions
function [N] = Common_Neighbors(vec1, vec2)
    N = sum(vec1 .* vec2);
end

function [N] = Jaccard_index(vec1, vec2)
    N = sum(vec1 .* vec2) / sum(ceil(vec1 + vec2)/2);
end

function [N] = Addamic_Addar(vec1, vec2)
    N = sum(1/log(sum(vec1 .* vec2)));
end