%% Analysing FOXP3 
% Intersection, applying new intesity filters, plot dynamics of
% selected upstream genes, plot biographs
% July 2017

%% FOXP3 PROBE SET 2 / 3    

GeneOfInterest = 'FOXP3';
IndexesOfProbessToUse = [2,3];
% CutoffGenesNumber = CutoffGenesNumberFOXP3;
FileNameGeneOfInterestDynamics = [GeneOfInterest '_' CellType 'Cell12_Probes23.mat'];

for IndexGeneName = 1:size(IndexesOfProbessToUse,2)
    
    %SelectedGeneName = GeneOfInterestProbesInCell12{IndexGeneName};
    VariableNameToLoad = [GeneOfInterest '_Probe' num2str(IndexesOfProbessToUse(IndexGeneName)) '_' CellType '12']; 
    
    % load file with resulting ranks for Cell 1, current probe
    FileNameToLoad = [GeneOfInterest '_Results_' CellType '_AllGenes_Cell1ProbeSet_', ...
        num2str(IndexesOfProbessToUse(IndexGeneName)),'.mat'];        
    load(FileNameToLoad, 'ranks')
    ranks2_Treg1 = ranks(1:CutoffGenesNumberFOXP3C1,:);
 
    % load file with resulting ranks for Cell 2, current probe
    FileNameToLoad = [GeneOfInterest '_Results_' CellType '_AllGenes_Cell2ProbeSet_', ...
        num2str(IndexesOfProbessToUse(IndexGeneName)),'.mat'];        
    load(FileNameToLoad, 'ranks')
    ranks2_Treg2 = ranks(1:CutoffGenesNumberFOXP3C2,:);

    
    % XXXXXX INTERSECT/UNION PROBES
    ProbeNumber = IndexesOfProbessToUse(IndexGeneName);

    % Calling script to perform intersection in Cell 1 and 2 and then apply intensity filters
    [Filtered_CutIntersection_ProbeN,GeneNamesCell12, ProbeIDsCell1, PublicIDsCell1,...
     ProbeIDsCell2, PublicIDsCell2] = IntersectSelectAnalyse_Ranking...
     (ranks2_Treg1, ranks2_Treg2,ThresholdAverageIntensityCell1,...
     ThresholdAverageIntensityCell2,CellType);

    if ProbeNumber == 2
       Filtered_CutIntersection_Probe_2 = Filtered_CutIntersection_ProbeN;
    end
    if ProbeNumber == 3
       Filtered_CutIntersection_Probe_3 = Filtered_CutIntersection_ProbeN;
    end
    
end
%--------------------------

   % Loading post-filter genes directly to save computation time for 14712 x 14472
   load('Filtered_CutIntersection_Probe_2_C1')   % 12815 genes
   load('Filtered_CutIntersection_Probe_2_C2')   % 12815 genes
   load('Filtered_CutIntersection_Probe_3_C1')   % 12327 genes
   load('Filtered_CutIntersection_Probe_3_C2')   % 12327 genes

   TcellDataTable = Tcell2DataTable;
   TableOfExtractedGeneNameAndIDs = TableOfExtractedGeneNameAndIDs_cell2;
   Cell = '2';
   SelectedGeneData = [FOXP3_Probe2_Treg12(2,:);FOXP3_Probe3_Treg12(2,:)];
   
 
% -------------------------------------------------------------------------  
% MAIN BIG LISTS TO REFER - WITHOUT FOXP3 AVG. FILTER
  % Genes regulating FOXP3 Probe-2 in Cell-1 (column 1:6) and Cell-2 (column 7:12)
  Filtered_CutIntersection_Probe_2 = Filtered_CutIntersection_Probe_2_C2;
  % Genes regulating FOXP3 Probe-3 in Cell-1 (column 1:6) and Cell-2 (column 7:12)
  Filtered_CutIntersection_Probe_3 = Filtered_CutIntersection_Probe_3_C2;
% -------------------------------------------------------------------------


  
%% =====================   APPROACH - 1   =====================
   
% Adding columns of 'Cell-1/2' 'Probe-2/3' in the end
% CELL-1
Filtered_CutInfo_Probe_2 = [Filtered_CutIntersection_Probe_2(:,1:6),...
                      cellstr(repmat('Probe-2',size(Filtered_CutIntersection_Probe_2,1),1)),...
                      cellstr(repmat('C1',size(Filtered_CutIntersection_Probe_2,1),1)),...
                      Filtered_CutIntersection_Probe_2(:,7:12),...
                      cellstr(repmat('Probe-2',size(Filtered_CutIntersection_Probe_2,1),1)),...
                      cellstr(repmat('C2',size(Filtered_CutIntersection_Probe_2,1),1))];
   
   
Filtered_CutInfo_Probe_3 = [Filtered_CutIntersection_Probe_3(:,1:6),...
                      cellstr(repmat('Probe-3',size(Filtered_CutIntersection_Probe_3,1),1)),...
                      cellstr(repmat('C1',size(Filtered_CutIntersection_Probe_3,1),1)),...
                      Filtered_CutIntersection_Probe_3(:,7:12),...
                      cellstr(repmat('Probe-3',size(Filtered_CutIntersection_Probe_3,1),1)),...
                      cellstr(repmat('C2',size(Filtered_CutIntersection_Probe_3,1),1))];
   
   
CombineBigLists = [Filtered_CutInfo_Probe_2; Filtered_CutInfo_Probe_3];
   
   
% To find Genes regulating FOXP3 probe-2 or 3 BUT, it must regulate in both - CELL 1 AND 2 
% i.e. DISCARD gene if ProbeID in Cell-1 does not matches with that in Cell-2
   count = 0;
    for ii = 1:size(CombineBigLists,1)
        count = count+1;
        SearchSameProbeID(:,count) = strcmp(CombineBigLists(ii,5),CombineBigLists(ii,5+8));
        LocOfSameProbeID = find(SearchSameProbeID' == 1);
    end  
    SameProbeID_C1C2 = CombineBigLists(LocOfSameProbeID,:);
    
    
Cell1_List = SameProbeID_C1C2(:,1:8);
Cell2_List = SameProbeID_C1C2(:,9:16);

% This list contains genes that regulate either probe2 or probe3 BUT they
% regulate in both cells, 1 and 2
CombineC1C2 = [Cell1_List; Cell2_List];



    

%-------------------------------------------   
%% Ranking by fitness score in CELL-1 and 2
% % % Now we rank 'CombineC1C2'
% % NewRank = [];
% %    for k=1:size(CombineC1C2,2)   
% %         cn = 1;
% %         FileToLoad = CombineC1C2;
% %         lineNumber = (1:size(FileToLoad,1))';
% %         locationOfGenes = cell2mat(FileToLoad(lineNumber,cn));
% %         fit_score = cell2mat(FileToLoad(lineNumber,cn+1));
% %         names = FileToLoad(lineNumber,cn+2);
% %         ActRep = cell2mat(FileToLoad(lineNumber,cn+3));
% %         ProbeID1a = FileToLoad(lineNumber,cn+4);
% %         RepPubID = FileToLoad(lineNumber,cn+5);
% %         [ranks_Cell, idx] = ranking_list_FOXP3(fit_score,locationOfGenes,names,ActRep,ProbeID1a,RepPubID);
% %         Ranks_Cell =  FileToLoad(idx,:);
% % 
% %         Cell_RankedProbes = [NewRank Ranks_Cell];   
% %    end
   
   
   
%--------------------------
%% PLOT DYNAMICS AND DIGRAPHS
    
    if PlotBioGraph == 1     
        Plot_4Figures
%         Plot_4Figures_another
    end  
   
   
   
   

%% =====================   APPROACH - 2   =====================
% % 
% %    Set_Tcell_1 = {Filtered_CutIntersection_Probe_2_C1 Filtered_CutIntersection_Probe_3_C1};
% %    Set_Tcell_2 = {Filtered_CutIntersection_Probe_2_C2 Filtered_CutIntersection_Probe_3_C2};
% % 
% % % Genes regulating PROBE-2 in CELL 1 AND 2 
% %    count = 0;
% %     for ii = 1:size(Filtered_CutIntersection_Probe_2_C2,1)
% %         count = count+1;
% %         SearchCommon_P2(:,count) = strcmp(Filtered_CutIntersection_Probe_2_C2(ii,5),Filtered_CutIntersection_Probe_2_C2(ii,5+6));
% %         LocOfCommon_P2 = find(SearchCommon_P2' == 1);
% %     end  
% %     SameProbeID_Probe2 = Filtered_CutIntersection_Probe_2_C2(LocOfCommon_P2,:);
% %     
% %     
% % % Genes regulating PROBE-3 in CELL 1 AND 2 
% %    count = 0;
% %     for jj = 1:size(Filtered_CutIntersection_Probe_3_C2,1)
% %         count = count+1;
% %         SearchCommon_P3(:,count) = strcmp(Filtered_CutIntersection_Probe_3_C2(jj,5),Filtered_CutIntersection_Probe_3_C2(jj,5+6));
% %         LocOfCommon_P3 = find(SearchCommon_P3' == 1);
% %     end  
% %     SameProbeID_Probe3 = Filtered_CutIntersection_Probe_3_C2(LocOfCommon_P3,:);    
% %     
% % % CELL-1
% %    CombineProbesC1 = [SameProbeID_Probe2(:,1:6); SameProbeID_Probe3(:,1:6)];
% %    CombineStringsC1 = [cellstr(repmat('Probe-2',size(SameProbeID_Probe2,1),1));...
% %                      cellstr(repmat('Probe-3',size(SameProbeID_Probe3,1),1))];
% %    Cell_1_list = [CombineProbesC1, CombineStringsC1, cellstr(repmat('C1',size(CombineProbesC1,1),1))];
% % % CELL-2
% %    CombineProbesC2 = [SameProbeID_Probe2(:,7:12); SameProbeID_Probe3(:,7:12)];
% %    CombineStringsC2 = [cellstr(repmat('Probe-2',size(SameProbeID_Probe2,1),1));...
% %                      cellstr(repmat('Probe-3',size(SameProbeID_Probe3,1),1))];
% %    Cell_2_list = [CombineProbesC2, CombineStringsC2, cellstr(repmat('C2',size(CombineProbesC2,1),1))];   
% %    
% %    
% %    Set_Tcells = {Cell_1_list, Cell_2_list};
% % 
% %    
% % %-------------------------------------------   
% % %% Ranking by fitness score in CELL-1 and 2
% %    for k=1:size(Set_Tcell_1,2)   
% %        names_Ranked = {'Cell1_Rank', 'Cell2_Rank'};
% % 
% %         cn = 1;
% %         FileToLoad = Set_Tcells{k};        % Selecting each Probe with changing 'k'
% %         lineNumber = (1:size(FileToLoad,1))';
% %         locationOfGenes = cell2mat(FileToLoad(lineNumber,cn));
% %         fit_score = cell2mat(FileToLoad(lineNumber,cn+1));
% %         names = FileToLoad(lineNumber,cn+2);
% %         ActRep = cell2mat(FileToLoad(lineNumber,cn+3));
% %         ProbeID1a = FileToLoad(lineNumber,cn+4);
% %         RepPubID = FileToLoad(lineNumber,cn+5);
% %         [ranks_Cell, idx] = ranking_list_FOXP3(fit_score,locationOfGenes,names,ActRep,ProbeID1a,RepPubID);
% %         Ranks_Cell =  FileToLoad(idx,:);
% % 
% %         Cell_Rank.(names_Ranked{k}) = Ranks_Cell;
% %         Cell_RankedProbes = struct2cell(Cell_Rank);   
% %    end
% % 
% %    CombinedCells = [Cell_RankedProbes{1};Cell_RankedProbes{2}];  
% %    
% %  %=====================
% %   
% %  % Now combine Cell1 ranked and Cell-2 ranked lists then rank them according to fitness score  
% %    CombRank = [];
% %    for kk=1:size(CombinedCells,2)   
% % 
% %         cn = 1;
% %         FileToLoad = CombinedCells;        % Selecting each Probe with changing 'k'
% %         lineNumber = (1:size(FileToLoad,1))';
% %         locationOfGenes = cell2mat(FileToLoad(lineNumber,cn));
% %         fit_score = cell2mat(FileToLoad(lineNumber,cn+1));
% %         names = FileToLoad(lineNumber,cn+2);
% %         ActRep = cell2mat(FileToLoad(lineNumber,cn+3));
% %         ProbeID1a = FileToLoad(lineNumber,cn+4);
% %         RepPubID = FileToLoad(lineNumber,cn+5);
% %         [ranks_Cell, idxComb] = ranking_list_FOXP3(fit_score,locationOfGenes,names,ActRep,ProbeID1a,RepPubID);
% %         Ranks_Comb =  FileToLoad(idxComb,:);
% % 
% %         Combined_Ranked = [CombRank Ranks_Comb];
% %    end
% %    
% %    %=====================
% %    
% % % WRITE FOR LOOP TO Search for same ProbeIDs in Cell-1 and Cell-2
% %    SearchGene = strcmp(Cell_RankedProbes{1}(:,5),Cell_RankedProbes{2}(:,5));
   


   
   
   
   
% % %--------------------------
% % %% PLOT DYNAMICS AND DIGRAPHS
% %     
% %     if PlotBioGraph == 1     
% %         PlotDynamicsBiographs_Ranking
% %     end
% %    
% %     if PlotDiGraph == 1
% %         plot_Digraph(Filtered_CutIntersection,GeneOfInterest,ProbeNumber,CellType);
% %     end    
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      END      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Ranking by fitness score in CELL-1
%    for k=1:size(Set_Tcell_1,2)    
%         names_Ranked_C1 = {'Probe2_Cell1_Rank', 'Probe3_Cell1_Rank'};
% 
%         cn = 1;
%         FileToLoad = Set_Tcell_1{k};        % Selecting each Probe with changing 'k'
%         lineNumber = (1:size(FileToLoad,1))';
%         locationOfGenes = cell2mat(FileToLoad(lineNumber,cn));
%         fit_score = cell2mat(FileToLoad(lineNumber,cn+1));
%         names = FileToLoad(lineNumber,cn+2);
%         ActRep = cell2mat(FileToLoad(lineNumber,cn+3));
%         ProbeID1a = FileToLoad(lineNumber,cn+4);
%         RepPubID = FileToLoad(lineNumber,cn+5);
%         [ranks_Cell1, idx_C1] = ranking_list_FOXP3(fit_score,locationOfGenes,names,ActRep,ProbeID1a,RepPubID);
%         Ranks_Cell1 =  FileToLoad(idx_C1,:);
% 
%         Cell1_Rank.(names_Ranked_C1{k}) = Ranks_Cell1;
%         Cell1_RankedProbes = struct2cell(Cell1_Rank);   
%    end


