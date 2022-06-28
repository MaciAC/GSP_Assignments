%% 1. Reading and cleaning up the data

T = readtable('LazegaLawyers/ELwork36.dat');

%% 2. Plot the graph
A = table2array(T);
[g,nodenums] = binaryImageGraph(A,4);
xcoor = g.Nodes.x;
ycoor = size(nodenums,2)-g.Nodes.y; % Flip to proper plot
figure(1);
plotImageGraph(g)

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

for i = selected_idx
    Adjacency_orig(edges(i,1), edges(i,2)) = Common_Neighbors(A(edges(i,1),:), A(edges(i,2),:));
    Adjacency_orig(edges(i,2), edges(i,1)) = Adjacency_orig(edges(i,1), edges(i,2));
end

%% 4. Probability of detection and false alarm
figure(1)

for thr = 1:10

Adjacency = Adjacency_orig;
Adjacency(Adjacency_orig < thr) = 0.0;


Detection = Adjacency .* A;
temp = zeros(N_nodes);
temp(Adjacency > 0) = 1;
False_Alarm = A - temp;
n_detected = 0;
n_expected = 0;
n_false_alarm = 0;
j=1;
A_roc = zeros(width(selected_idx),10);
Adjacency_roc = zeros(width(selected_idx),10);

for i = selected_idx
    A_roc(j,thr) =  A(edges(i,1), edges(i,2));
    Adjacency_roc(j,thr) = Adjacency(edges(i,1), edges(i,2));
    if A(edges(i,1), edges(i,2)) > 0
        n_expected = n_expected + 1;
    end
    if Detection(edges(i,1), edges(i,2)) > 0
        n_detected = n_detected + 1; 
    end
    if False_Alarm(edges(i,1), edges(i,2)) == -1
        n_false_alarm = n_false_alarm + 1;
    end
    j = j +1;
end

P_detection = n_detected / n_expected

P_false_alarm = n_false_alarm / n_expected

plotroc(reshape(A_roc,1,[]), reshape(Adjacency_roc,1,[]))
hold on
end

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


