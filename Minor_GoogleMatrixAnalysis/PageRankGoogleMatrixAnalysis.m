
clear; close all; clc;

FitnessThrashold = 70.0;
DampingFactor = 0.85;
load('Results_Teff_MIGHT_Cell2.mat', 'fitness');

%%%%%%%% TRANSFORM FITNESS TO AN ADJACENCY MATRIX %%%%%%%
figure;
image(fitness)
title('Fitness Matrix'); xlabel('Genes'); ylabel('Genes');

colorbar;

MatrixDimension = size(fitness,1);
NetworkTopology = cell(MatrixDimension,MatrixDimension);
for i = 1:MatrixDimension
    for j = 1:MatrixDimension
        if fitness(i,j) >= FitnessThrashold 
            NetworkTopology{i,j} = 1;
        elseif fitness(i,j) < FitnessThrashold 
            NetworkTopology{i,j} = 0;
        else
            disp('Error!')
        end
    end
end
AdjacencyMatrix = cell2mat(NetworkTopology); % fitness/100.;%
figure;
spy(AdjacencyMatrix);
title('Adjacency Matrix');  xlabel('Genes'); ylabel('Genes');

G = digraph(AdjacencyMatrix);

%%%%%%%% NOW COMPUTE GOOGLE MATRIX AND PAGE RANK %%%%%%%
format long
pr = centrality(G,'pagerank','FollowProbability',DampingFactor);
G.Nodes.PageRank = pr;
G.Nodes.InDegree = indegree(G);
G.Nodes.OutDegree = outdegree(G);

G.Nodes

figure;
plot(G,'NodeLabel',{},'NodeCData',G.Nodes.PageRank,'Layout','force');
title('Directed Graph, Nodes Colored by Page Rank');
colorbar;

figure;
H = subgraph(G,find(G.Nodes.PageRank > 0.005));
plot(H,'NodeLabel',{},'NodeCData',H.Nodes.PageRank,'Layout','force');
title('Directed Graph, only nodes with higher Page Rank, Nodes Colored by Page Rank');
colorbar;

SortedFitnessLong = sortrows([G.Nodes.PageRank,fitness']);
SortedFitness = SortedFitnessLong(:,2:MatrixDimension+1)';

figure;
image(SortedFitness)
title('Fitness Matrix, columns sorted by PageRank'); xlabel('Genes'); ylabel('Genes');
colorbar;

%%%%%%%% NOW SAME BUT WITH TRANSPOSE OF ADJACENCY MATRIX TO GET THE CHEI RANK %%%%%%%

AdjacencyMatrixOfComplement = AdjacencyMatrix';

GI = digraph(AdjacencyMatrixOfComplement);

prI = centrality(GI,'pagerank','FollowProbability',DampingFactor);
GI.Nodes.PageRank = prI;
GI.Nodes.InDegree = indegree(GI);
GI.Nodes.OutDegree = outdegree(GI);
GI.Nodes

figure;
plot(GI,'NodeLabel',{},'NodeCData',GI.Nodes.PageRank,'Layout','force');
title('Directed Graph, Nodes Colored by Chei Rank');
colorbar;

figure;
HI = subgraph(GI,find(GI.Nodes.PageRank > 0.01));
plot(HI,'NodeLabel',{},'NodeCData',HI.Nodes.PageRank,'Layout','force');
title('Directed Graph, only nodes with higher Chei Rank, Nodes Colored by Chei Rank');
colorbar;

%%%%%%%%% Plotting into PageRank-CheiRank plane %%%%%%%%%

xy = [G.Nodes.PageRank GI.Nodes.PageRank];
figure;
plot(G,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',5);
title('Directed Graph'); xlabel('PageRank'); ylabel('CheiRank');

rnkPR = floor(tiedrank(G.Nodes.PageRank));
x = (MatrixDimension+1)*ones(MatrixDimension,1)-rnkPR;

rnkCR = floor(tiedrank(GI.Nodes.PageRank));
y = (MatrixDimension+1)*ones(MatrixDimension,1)-rnkCR;

xyInd = [x y];

figure;
plot(G,'XData',xyInd(:,1),'YData',xyInd(:,2),'MarkerSize',5);
title('Directed Graph'); xlabel('Rank in PageRank ranking'); ylabel('Rank in CheiRank ranking');

figure;
scatter(xyInd(:,1), xyInd(:,2))
title('Nodes'); xlabel('Rank in PageRank ranking'); ylabel('Rank in CheiRank ranking');

figure;
hist3([xyInd(:,1), xyInd(:,2)],[12 12]);
colorbar;
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('Density of Nodes'); xlabel('Rank in PageRank ranking'); ylabel('Rank in CheiRank ranking');
grid on;
view(2);

%%%%%%%%% Plotting in LogLog scale into PageRank-CheiRank plane %%%%%%%%%
figure;
scatter(xyInd(:,1), xyInd(:,2))
title('Nodes'); xlabel('Log(Rank in PageRank ranking)'); ylabel('Log(Rank in CheiRank ranking)');
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');

figure;
hist3([xyInd(:,1), xyInd(:,2)],[12 12]);
colorbar;
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
title('Density of Nodes'); xlabel('Log(Rank in PageRank ranking)'); ylabel('Log(Rank in CheiRank ranking)');
grid on;
view(2);

%%%%%%%%% Compute Page Rank - Chei Rank Correlation %%%%%%%%%% 
N = MatrixDimension
k = N * sum((G.Nodes.PageRank .* GI.Nodes.PageRank)) - 1


