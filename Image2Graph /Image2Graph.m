clear all;close all
run('/Users/alba/Documents/ETSETB/GRAPH SP/gspbox/gsp_start.m')
run('/Users/alba/Documents/ETSETB/GRAPH SP/unlocbox/init_unlocbox.m')
%run('/Users/alba/Documents/ETSETB/GRAPH SP/unlocbox/init_unlocbox.m')

% Look for an image of 50x50 pixels, increase contrast to achieve a
% connected graph

bw=imread('snake.png');
bw =sign(imcomplement(bw));
bw=bw(:,:,1);
% [i,j]=find(bw==0);
% for ind=1:length(i)
% bw(i(ind),j(ind))=1;
% end;
[i,j]=find(bw>1);
for ind=1:length(i)
bw(i(ind),j(ind))=1;
end;
[g,nodenums] = binaryImageGraph(bw,4);
xcoor = g.Nodes.x;
ycoor = size(nodenums,2)-g.Nodes.y; % Flip to proper plot
figure(2);
plotImageGraph(g)

