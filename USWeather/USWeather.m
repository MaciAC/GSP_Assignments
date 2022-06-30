run('/Users/maciaac/Documents/UPC/2n Semestre/GSP/GSP_Assignments/gspbox/gsp_start.m');
%% Reading and visualization of the data

A = load("US Weather Station Data/adjacency.mat").A;
position = load("US Weather Station Data/position.mat").position_matrix;
T_Farenheit = load("US Weather Station Data/temp.mat").Y;
T_Celsius = (T_Farenheit - 32) * 5/9;

%% Generation of an undirected graph
figure(1)
mesh(T_Celsius);

figure(2)
G = gsp_graph(A,position);
gsp_plot_graph(G);

hawaii_idxs = find(position(:,1)  < 24);

A(hawaii_idxs,:) = [];
A(:,hawaii_idxs) = [];
position(hawaii_idxs,:) = []; 
T_Celsius(hawaii_idxs,:) = [];

figure(3)
G = gsp_graph(A,position);
gsp_plot_graph(G);

figure(4)
A = A + A';
G = gsp_graph(A,position);
gsp_plot_graph(G);

%%
