%% INTERSECTING GENES IN THE RANKED LIST OF Treg-1 and Treg-2

clc;
clear;
close all;

%---------------
N = 500;
SelectCandidates = 21;
FromCell_1 = 0;
FromCell_2 = 1;

PlotFigure = 1;
%---------------
        
load('FILTERED_ExtractedCutData_TregALLgen.mat')
load('FOXP3_TregCell12_Probes23.mat')
load('Cell1_List.mat')
load('Cell2_List.mat')
load('CodingNonCodingList.mat');
% load('Final_P2P3_rank_Treg1.mat')
% load('Final_P2P3_rank_Treg2.mat')

CutOff_C1 = size(Cell1_List,1);
CutOff_C2 = size(Cell2_List,1);

CodingNonCodSHORTlist = [];
time_points = 0:20:360;


%---------------
%% Ranking by fitness score in CELL-1 and 2
% Ranking cell 1
NewRank = [];
%    for k=1:size(Cell1_List,1)   
        cn = 1;
        FileToLoad = Cell1_List;
        lineNumber = (1:size(FileToLoad,1))';
        locationOfGenes = cell2mat(FileToLoad(lineNumber,cn));
        fit_score = cell2mat(FileToLoad(lineNumber,cn+1));
        names = FileToLoad(lineNumber,cn+2);
        ActRep = cell2mat(FileToLoad(lineNumber,cn+3));
        ProbeID1a = FileToLoad(lineNumber,cn+4);
        RepPubID = FileToLoad(lineNumber,cn+5);
        [ranks_Cell, idx] = ranking_list_FOXP3(fit_score,locationOfGenes,names,ActRep,ProbeID1a,RepPubID);
        Ranks_Cell =  FileToLoad(idx,:);

        Cell_RankedProbes = [NewRank Ranks_Cell];   
%    end
   
   Cell1_List_Ranked = Cell_RankedProbes;

%---------------
    ColumnNumberForIntersect = 3; % 3 means use gene name to intersect
    IntersectionTreg1Treg2 = intersectRanks_V2(Cell1_List_Ranked(1:N,:),...
                             Cell2_List(1:N,:),...
                             ColumnNumberForIntersect);

    % Keep lines with gene name = "---" if and only if they have same
    % public id for both cell 1 and cell 2
    FindDashGenes = strcmp('---',IntersectionTreg1Treg2(:,3));
    LocationsOfNormalGenes = find(FindDashGenes == 0);
    LocationsOfDashGenes = find(FindDashGenes == 1);
    
    PublicIDsC1 = IntersectionTreg1Treg2(:,6);
    PublicIDsC2 = IntersectionTreg1Treg2(:,6+8);
    ComparePublicIDsC12 = strcmp(PublicIDsC1,PublicIDsC2);
    IsNotDashDash = not(strcmp(IntersectionTreg1Treg2(:,3),'---'));
    MyNewVector = (IsNotDashDash == 1 | ComparePublicIDsC12 == 1);
    findLocationOfCompare = find(MyNewVector);
    LinesToBeKept = IntersectionTreg1Treg2(findLocationOfCompare,:);

    % We select only those lines for which in both cells/probes upstream
    % regulating genes are ALWAYS REPRESSING or ALWAYS ACTIVATING
    temp_Probe2 = find(cell2mat(LinesToBeKept(:,4)) ...
        ./ cell2mat(LinesToBeKept(:,4+6)) > 0);
    CutIntersectionMatrix_Treg1Treg2 = LinesToBeKept(temp_Probe2,:);

    % Remove duplicates
    [Filtered_Cell1_wd] = remove_duplicates(CutIntersectionMatrix_Treg1Treg2(:,1:8)); 
    [Filtered_Cell2_wd] = remove_duplicates(CutIntersectionMatrix_Treg1Treg2(:,9:16)); 
    
    % Ranking by fitness score
    Final_rank_C1 = [];
    cn = 1;
    FileToLoad = Filtered_Cell1_wd(:,:);
    lineNumber = (1:size(FileToLoad,1))';
    locationOfGenes = cell2mat(FileToLoad(lineNumber,cn));
    fit_score = cell2mat(FileToLoad(lineNumber,cn+1));
    names = FileToLoad(lineNumber,cn+2);
    ActRep = cell2mat(FileToLoad(lineNumber,cn+3));
    ProbeID1a = FileToLoad(lineNumber,cn+4);
    RepPubID = FileToLoad(lineNumber,cn+5);
    ProbeFOXP3 = FileToLoad(lineNumber,cn+6);
    MainCell = FileToLoad(lineNumber,cn+7);
    [ranks_Cell1, idx] = ranking_list_FOXP3(fit_score,locationOfGenes,...
                         names,ActRep,ProbeID1a,RepPubID,ProbeFOXP3,MainCell);
    Ranks_Cell1 =  FileToLoad(idx,:);  

    Filtered_rank_C1 = [Final_rank_C1,Ranks_Cell1];
    Filtered_Treg1_rank = Filtered_rank_C1;
    
    %----------------
    Final_rank_C2 = [];
    cn = 1;
    FileToLoad = Filtered_Cell2_wd(:,:);
    lineNumber = (1:size(FileToLoad,1))';
    locationOfGenes = cell2mat(FileToLoad(lineNumber,cn));
    fit_score = cell2mat(FileToLoad(lineNumber,cn+1));
    names = FileToLoad(lineNumber,cn+2);
    ActRep = cell2mat(FileToLoad(lineNumber,cn+3));
    ProbeID1a = FileToLoad(lineNumber,cn+4);
    RepPubID = FileToLoad(lineNumber,cn+5);
    ProbeFOXP3 = FileToLoad(lineNumber,cn+6);
    MainCell = FileToLoad(lineNumber,cn+7);
    [ranks_Cell2, idx] = ranking_list_FOXP3(fit_score,locationOfGenes,...
                         names,ActRep,ProbeID1a,RepPubID,ProbeFOXP3,MainCell);
    Ranks_Cell2 =  FileToLoad(idx,:);  

    Filtered_rank_C2 = [Final_rank_C2,Ranks_Cell2];
    Filtered_Treg2_rank = Filtered_rank_C2;
    
    GenesToCheck(1:SelectCandidates,9) = num2cell(zeros(size(GenesToCheck,1),1));
    
% Checking for the same ProbeIDs   
    for ii = 1:size(Filtered_Treg1_rank,1)
        if cell2mat(Filtered_Treg1_rank(ii,9)) == 0
            
        GeneSearchInC1C2 = GenesToCheck(ii,3);
        SearchCommon_C1 = strcmp(Filtered_Treg1_rank(ii,5),Filtered_Treg2_rank(:,5));
        LocOfCommon_C1 = find(SearchCommon_C1 == 1);
        V_Cell1{ii} = LocOfCommon_C1; 
        InfoOnPlot_C1 = [Filtered_Treg1_rank(ii,:);...
                         Filtered_Treg2_rank(LocOfCommon_C1,:)];
        LocToPlot_C1 = cell2mat(Filtered_Treg1_rank(ii,1));   % Location in the main data file
        DataToPlot_C1 = Tcell1DataTable(LocToPlot_C1,:);
        end
    end
    
        
% %     for ii = 1:size(Filtered_Treg2_rank,1)
% %         SearchCommon_C2 = strcmp(Filtered_Treg2_rank(ii,5),Filtered_Treg1_rank(:,5));
% %         LocOfCommon_C2 = find(SearchCommon_C2 == 1);
% %         V_Cell2{ii} = LocOfCommon_C2;
% %         
% %         InfoOnPlot_C2 = Filtered_Treg2_rank(ii,:);
% %         LocToPlot_C2 = cell2mat(Filtered_Treg2_rank(ii,1));   % Location in the main data file
% %         DataToPlot_C2 = Tcell2DataTable(LocToPlot_C2,:);
% %     end
    
    % Collecting info and data-to-plot from the locations where the genes is found
%     InfoOnPlot_C1 = Filtered_Treg2_rank(LocOfCommon_C1,:)
%     LocToPlot_C1 = cell2mat(Filtered_Treg2_rank(LocOfCommon_C1,1));   % Location in the main data file
%     DataToPlot_C1 = Tcell1DataTable(LocToPlot_C1,:);
    
%     InfoOnPlot_C2 = Filtered_Treg1_rank(LocOfCommon_C2,:);
%     LocToPlot_C2 = cell2mat(Filtered_Treg1_rank(LocOfCommon_C2,1));   % Location in the main data file
%     DataToPlot_C2 = Tcell1DataTable(LocToPlot_C2,:);
    
    
%---------------





