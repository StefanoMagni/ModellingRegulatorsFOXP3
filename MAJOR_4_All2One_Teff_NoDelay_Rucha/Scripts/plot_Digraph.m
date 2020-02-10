 %% DIGRAPH
 
 
 function plot_Digraph(Filtered_CutIntersection_ProbeN,GeneOfInterest,ProbeNumber,CellType)
 
 
for Tregcell = 1:2
     if Tregcell == 1          
        colNo = 0;
        cellNumber = 1;  % no. of cells
     else
        colNo = 6;       % jump by 6 columns to reach column related to cell 2
        cellNumber = 2;  % no. of cells
     end     
    
    % TO REMANE REPEATED GENE NAMES AS, FOR EG:, ERC, ERC-2, ERC-3..AND SO ON
    % Here, TargetGene is one of the 4 genes of interest
    TargetGene = [num2str(GeneOfInterest) '-Probe' num2str(ProbeNumber)];    
    speciesC1 = {TargetGene,Filtered_CutIntersection_ProbeN{:,3+colNo}}';
    
    for i = 1:1:numel(speciesC1(:,:))
        k = 1;
        HasTwins = 0;
        TentativeNodeName = speciesC1{i,1};
        for j = 1:(i-1)
            if strcmp(TentativeNodeName,speciesC1{j})
                k = k + 1;
                HasTwins = 1;
            end
        end
        if HasTwins
           TentativeNodeName = [sprintf(TentativeNodeName), '-', num2str(k)];
        else
           TentativeNodeName = [sprintf(TentativeNodeName)];
        end
        mRNAs_names{i} = TentativeNodeName;
    end

    clear species
    speciesC1 = mRNAs_names;

    % PLOT DIGRAPH HERE
    sources = 1+1:size(Filtered_CutIntersection_ProbeN,1)+1;      
    target = repelem(1,size(Filtered_CutIntersection_ProbeN,1));
    names = {speciesC1{:}};
    weight = round(cell2mat(Filtered_CutIntersection_ProbeN(:,2+colNo)))';
    D = digraph(sources,target,weight,names);

    % PLOT DIGRAPH
    figure();
%     Digraph = plot(D,'ArrowSize',12);                      % for normal digraph
    Digraph = plot(D,'Layout','circle','ArrowSize',18);      % for circular digraph
    title([TargetGene '- ' CellType num2str(cellNumber)])
    % To enlarge figure to full screen.
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

    % Edit Digraph properties
    Digraph.MarkerSize = 2.5 ;
    Digraph.LineWidth = 1.2;
    Digraph.EdgeLabel = table2array(D.Edges(:,2)); % num2str(weight);

    full(adjacency(D));   % shows all2one matrix
    outdegree(D);         % shows number of edges for the node

    % Find the location of genes that are activating and repressing to color the edges
    DetectActivation = (cell2mat(Filtered_CutIntersection_ProbeN(:,4+colNo)) >= 0);
    Activation = find(DetectActivation == 1);    % Location of genes ACTIVATING FOXP3
    Repression = find(DetectActivation == 0);    % Location of genes REPRESSING FOXP3

    % HIGHLIGHT nodes and edges
    highlight(Digraph,Activation+1,'NodeColor',[0.6,0.9,0.6])  % highlighting nodes that are activating
    highlight(Digraph,Repression+1,'NodeColor',[1 0.4 0.6])  % highlighting nodes that are repressing

    highlight(Digraph,sources(Activation),target(Activation),...
             'EdgeColor','g','LineWidth',1.5)
    highlight(Digraph,sources(Repression),target(Repression),...
             'EdgeColor','r','LineWidth',1.5)
    highlight(Digraph,names)

    % EDIT the node label size
%     labelColors = parula(length(names));
    labelSizes  = repelem(8,size([names weight],2));
    Digraph.NodeLabel = {};
    % Custom labels
    hold on
    for i=1:length(names)
    text(Digraph.XData(i), Digraph.YData(i), names(i),...
        'Color', 'k' , 'FontSize', labelSizes(i));
    end
    hold off
    
    % SAVE DIGRAPHS HERE
    temp=[ num2str(GeneOfInterest) '_Probe' num2str(ProbeNumber) '_',...
           CellType num2str(cellNumber) '_DIGRAPH','.eps']; saveas(gcf,temp,'epsc');
    
    
end

    
    
end




    
    