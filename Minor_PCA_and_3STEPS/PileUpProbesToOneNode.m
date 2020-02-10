
%%%%%%%%% MERGE DIFFERENT PROBES IN ONE NODE %%%%%%%%%%

function [CUTmRNA_names, CutCutNewMODfitness, CutCututdcgains] = PileUpProbesToOneNode( ...
    TableOfExtractedGeneNameAndIDs_cellX, mRNA_namesFromLoadedFile, fitness, dcgains)

    MyTable = TableOfExtractedGeneNameAndIDs_cellX;
    mRNA_names = {MyTable{:,4}};
    mRNAs_names = {MyTable{1,4}};
    MyTable{1,6} = 1;
    StartingIndex = 1;
    for i = 2:1:numel(MyTable(:,4))
         k = 1;
         HasTwins = 0;
         TentativeNodeName = MyTable{i,4};
         for j = 1:(i-1)
             if strcmp(TentativeNodeName,mRNA_names{j})
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
         MyTable{i,6} = k;
         if not(HasTwins)
             StartingIndex = i;
             MyTable{i,7} = [StartingIndex:i];
         elseif HasTwins
         	MyTable{i,7} = [StartingIndex:i];
         end
    end
    
    MyStr = string(MyTable(:,4));
    [C,ia,ic] = unique(MyStr,'legacy');
    MyCutTable = MyTable(ia,:);
    
    %     CellIndex = 2;  % HEY!!!
    %     if CellIndex == 1                    
    %         ResultsFileNameToBeLoaded = ['3STEPS_III_Results_Treg_KNOWN_Cell1.mat'];
    %     elseif CellIndex == 2
    %         ResultsFileNameToBeLoaded = ['3STEPS_III_Results_Treg_KNOWN_Cell2.mat'];
    %     end
    %     load(ResultsFileNameToBeLoaded);

    % mRNA_names; 
    CUTmRNA_names = mRNA_namesFromLoadedFile(ia);

    % fitness;
    MODfitness = fitness;
    %%% Pile up LINES for same prob
    for i = 1:numel(MyCutTable(:,1))
        temp = MyCutTable(i,7);
        IndexesForCurrentProbe = temp{1};
        LinesForOneProbe = fitness(IndexesForCurrentProbe,:);
        if size(LinesForOneProbe,1)==1
            LineOfMaxima = LinesForOneProbe;
        elseif not(size(LinesForOneProbe,1)==1)
        	LineOfMaxima = max(LinesForOneProbe);
        end
        RepeatedLine = 0;
        for index = IndexesForCurrentProbe
            if not(RepeatedLine)
                MODfitness(index,:) = LineOfMaxima;
                RepeatedLine = 1;
            elseif RepeatedLine
                MODfitness(index,:) = -1*ones(1,numel(LineOfMaxima));
            end
        end
    end
    
    %%% Pile up COLUMNS for same prob
    TransposeOfMODfitness = MODfitness';
    NewMODfitness = TransposeOfMODfitness;
    %%% Pile up LINES for same prob
    for i = 1:numel(MyCutTable(:,1))
        temp = MyCutTable(i,7);
        IndexesForCurrentProbe = temp{1};
        LinesForOneProbe = TransposeOfMODfitness(IndexesForCurrentProbe,:);
        if size(LinesForOneProbe,1)==1
            LineOfMaxima = LinesForOneProbe;
        elseif not(size(LinesForOneProbe,1)==1)
        	LineOfMaxima = max(LinesForOneProbe);
        end
        RepeatedLine = 0;
        for index = IndexesForCurrentProbe
            if not(RepeatedLine)
                NewMODfitness(index,:) = LineOfMaxima;
                RepeatedLine = 1;
            elseif RepeatedLine
                NewMODfitness(index,:) = -1*ones(1,numel(LineOfMaxima));
            end
        end
    end

    Transposed = NewMODfitness';
    %%% Cut out useless lines
    CutNewMODfitness = Transposed(not(Transposed(:,1)==-1),:);
    %%% Cut out useless columns
    CutCutNewMODfitness = CutNewMODfitness(:,not(CutNewMODfitness(1,:)==-1));
    
    %%%%% NOW SAME ON dcgains!!! %%%%%
    %%% Cut out useless lines
    Cutdcgains = dcgains(not(Transposed(:,1)==-1),:);
    %%% Cut out useless columns
    CutCututdcgains = Cutdcgains(:,not(CutNewMODfitness(1,:)==-1));

    %%%%%%%%% END of MERGE DIFFERENT PROBES IN ONE NODE %%%%%%%%%%
end