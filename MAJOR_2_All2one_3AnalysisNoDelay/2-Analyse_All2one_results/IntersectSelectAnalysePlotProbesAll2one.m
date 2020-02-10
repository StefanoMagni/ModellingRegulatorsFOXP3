
function [Filtered_CutIntersection_ProbeN,GeneNamesCell12, ProbeIDsCell1, PublicIDsCell1,...
         ProbeIDsCell2, PublicIDsCell2] = IntersectSelectAnalysePlotProbesAll2one...
         (ranks2_Treg1, ranks2_Treg2,ThresholdAverageIntensityCell1,...
         ThresholdAverageIntensityCell2)
    
    load('FILTERED_ExtractedCutData_TregALLgen.mat')

    % Extracting Probe ID for the ranks above
    indexCell1 = cell2mat(ranks2_Treg1(:,1));
    indexCell2 = cell2mat(ranks2_Treg2(:,1));
    ExtractProbeIDCell1 = TableOfExtractedGeneNameAndIDs_cell1(indexCell1,1:2);
    ExtractProbeIDCell2 = TableOfExtractedGeneNameAndIDs_cell2(indexCell2,1:2);

    % Ranks with Probe ID and Public ID added in the last column no.5
    ranks2_Treg1 = [ranks2_Treg1(:,1:4) ExtractProbeIDCell1];
    ranks2_Treg2 = [ranks2_Treg2(:,1:4) ExtractProbeIDCell2];

    % Intersecting the genes from both the cells using PUBLIC ID
%     ColumnNumberForIntersect = 5; % 5 means use probe IDs to intersect
%     ColumnNumberForIntersect = 6; % 6 means use public IDs to intersect
    ColumnNumberForIntersect = 3; % 3 means use gene name to intersect
    IntersectionTreg1Treg2_Probe2 = intersectRanks(ranks2_Treg1,ranks2_Treg2,...
    ColumnNumberForIntersect);

    % Keep lines with gene name = "---" if and only if they have same
    % public id for both cell 1 and cell 2
%     FindDashGenes = strcmp('---',IntersectionTreg1Treg2_Probe2(:,3));
%     LocationsOfNormalGenes = find(FindDashGenes == 0);
%     LocationsOfDashGenes = find(FindDashGenes == 1);
%     
%     PublicIDsC1 = IntersectionTreg1Treg2_Probe2(LocationsOfDashGenes,6);
%     PublicIDsC2 = IntersectionTreg1Treg2_Probe2(LocationsOfDashGenes,6+6);
%     ComparePunlicIDsC12 = strcmp(PublicIDsC1,PublicIDsC2);
%     findLocationOfCompare = find(ComparePunlicIDsC12 == 1);
%     DashDashDashOnly = IntersectionTreg1Treg2_Probe2(findLocationOfCompare,:);

    FindDashGenes = strcmp('---',IntersectionTreg1Treg2_Probe2(:,3));
    LocationsOfNormalGenes = find(FindDashGenes == 0);
    LocationsOfDashGenes = find(FindDashGenes == 1);
    
    PublicIDsC1 = IntersectionTreg1Treg2_Probe2(:,6);
    PublicIDsC2 = IntersectionTreg1Treg2_Probe2(:,6+6);
    ComparePunlicIDsC12 = strcmp(PublicIDsC1,PublicIDsC2);
    IsNotDashDash = not(strcmp(IntersectionTreg1Treg2_Probe2(:,3),'---'));
    MyNewVector = (IsNotDashDash == 1 | ComparePunlicIDsC12 == 1);
    findLocationOfCompare = find(MyNewVector);
    LinesToBeKept = IntersectionTreg1Treg2_Probe2(findLocationOfCompare,:);

    % We select only those lines for which in both cells/probes upstream
    % regulating genes are ALWAYS REPRESSING or ALWAYS ACTIVATING
    temp_Probe2 = find(cell2mat(LinesToBeKept(:,4)) ...
        ./ cell2mat(LinesToBeKept(:,4+6)) > 0);
    CutIntersectionMatrix_Treg1Treg2_Probe2 = LinesToBeKept(temp_Probe2,:);

    % Extracting data of the filtered genes from the main file
    LocationsOfPublicIDCell1 = cell2mat(CutIntersectionMatrix_Treg1Treg2_Probe2(:,1));
    GenesToBeFilteredCell1 = Tcell1DataTable(LocationsOfPublicIDCell1,:);
    LocationsOfPublicIDCell2 = cell2mat(CutIntersectionMatrix_Treg1Treg2_Probe2(:,1+6));
    GenesToBeFilteredCell2 = Tcell2DataTable(LocationsOfPublicIDCell2,:);

    filter_1 = ThresholdAverageIntensityCell1;  % set filter threshold value
    filter_2 = ThresholdAverageIntensityCell2;

    avgCell1 = mean(GenesToBeFilteredCell1')';  % compute average for each gene
    avgCell2 = mean(GenesToBeFilteredCell2')';

    % Here results are: 0=below threshold and 1=above threshold
    avgFilteredGenesCell1 = find(avgCell1 >= filter_1);
    avgFilteredGenesCell2 = find(avgCell2 >= filter_2);

    % Provides common genes passing the filter in both the cells
    intersectFilteredGenes = intersect(avgFilteredGenesCell1,avgFilteredGenesCell2);     
    Filtered_CutIntersection_ProbeN = ...
                CutIntersectionMatrix_Treg1Treg2_Probe2(intersectFilteredGenes,:);

    % Create an array of strings to display in figures
    dispProbeID = repmat('ProbeID=',size(Filtered_CutIntersection_ProbeN,1),1);
    dispPublicID = repmat('PublicID=',size(Filtered_CutIntersection_ProbeN,1),1);
    dispGeneName = repmat('Gene Name=',size(Filtered_CutIntersection_ProbeN,1),1);

    % To display Probe ID - Public ID - Gene Name in figures           
    GeneNamesCell12 = [dispGeneName char(Filtered_CutIntersection_ProbeN(:,3))];

    ProbeIDsCell1 = [dispProbeID char(Filtered_CutIntersection_ProbeN(:,5))];
    PublicIDsCell1 = [dispPublicID char(Filtered_CutIntersection_ProbeN(:,6))];

    ProbeIDsCell2 = [dispProbeID char(Filtered_CutIntersection_ProbeN(:,5+6))];
    PublicIDsCell2 = [dispPublicID char(Filtered_CutIntersection_ProbeN(:,6+6))];

end


 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    END    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end