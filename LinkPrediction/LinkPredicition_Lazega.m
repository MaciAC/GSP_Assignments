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
X = randi(N_edges) - 1;
Y = randi(N_nodes);
