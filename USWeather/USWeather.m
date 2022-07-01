run('/Users/maciaac/Documents/UPC/2n Semestre/GSP/GSP_Assignments/gspbox/gsp_start.m');
%% Reading and visualization of the data

A = load("US Weather Station Data/adjacency.mat").A;
pos = load("US Weather Station Data/position.mat").position_matrix;
T_Farenheit = load("US Weather Station Data/temp.mat").Y;
T_Celsius = (T_Farenheit - 32) * 5/9;

position = pos;
position(:,1) = pos(:,2);
position(:,2) = pos(:,1);


figure(1)
mesh(T_Celsius);

%% Generation of an undirected graph
figure(2)
G = gsp_graph(A,position);

gsp_plot_graph(G);

hawaii_idxs = find(position(:,2)  < 24);

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

%% Low pass filter
% shift operator (laplacian)
D = zeros(height(A));
for i = 1:height(A)
    D(i,i) = sum(A(i,:));
end

L = D - A;

[Vec, Val] = eig(L);
x = T_Celsius(:,1);
x_ = conj(Vec)'*x;

h_lowpass = zeros(height(x_));
for i = 1:height(x_)
    h_lowpass = 1 / (1 + Val(i,i));
end
y_ = x_ .* h_lowpass;

y_l = conj(Vec)'*y_;


%% High pass filter
x = T_Celsius(:,1);
x_ = conj(Vec)'*x;

h_highpass = zeros(height(x_));
max_eigenval = max(max(Val));
for i = 1:height(x_)
    h_highpass = 1 / (1 - (Val(i,i)- max_eigenval));
    
end
y_ = x_ .* h_highpass;

y_h = conj(Vec)'*y_;

figure(10)
mesh([x,y_l])
figure(11)
mesh([x,y_h])

%% different seasons
days = [1, 91, 181, 271];
max_eigenval = max(max(Val));
for d = days
    x = T_Celsius(:,d);
    x_ = conj(Vec)'*x;
    
    h_lowpass = zeros(height(x_),1);
    h_highpass = zeros(height(x_),1);
    for i = 1:height(x_)
        h_highpass(i) = 1 / (1 - (Val(i,i) - max_eigenval));
        h_lowpass(i) = 1 / (1 + Val(i,i));
    end
    y_high = x_ .* h_highpass;
    y_low = x_ .* h_lowpass;
    
    y_h = Vec*y_high;
    y_l = Vec*y_low;
    figure(d);
    G = gsp_graph(A,position);
    %param.climits=[0 1];
    subplot(211)
    gsp_plot_signal(G,y_l);
    subplot(212)
    gsp_plot_signal(G,y_h);
end