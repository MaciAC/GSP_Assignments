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
X = randi(N_nodes, N_edges/5, 2);
for i = 1:height(X)
    if X(i,1) == X(i,2)
        if X(i,2) == N_nodes
            X(i,2) = 1;
        else
            X(i,2) = X(i,2) + 1;
        end
    end
    Adjacency_orig(X(i,1), X(i,2)) = Common_Neighbors(A(X(i,1),:), A(X(i,2),:));
    Adjacency_orig(X(i,2), X(i,1)) = Adjacency_orig(X(i,1), X(i,2));
end


%% 4. Probability of detection and false alarm
for thr = 1:10
figure(thr)

Adjacency = zeros(N_nodes);
Adjacency(Adjacency_orig < thr) = 0;
Adjacency(Adjacency_orig > 0) = 1;

Detection = Adjacency .* A;
False_Alarm = A - Adjacency;

n_detected = sum(Detection(Detection > 0)) / 2;
P_detection = n_detected / (sum(A(A > 0))/2);

n_false_alarm = sum(False_Alarm(False_Alarm == -1)) / 2;
P_false_alarm = n_false_alarm / (sum(A(A > 0))/2);

plotroc(reshape(A,1,[]), reshape(Adjacency,1,[]))
end

%% functions
function [N] = Common_Neighbors(vec1, vec2)
    N = sum(vec1 .* vec2);
end


