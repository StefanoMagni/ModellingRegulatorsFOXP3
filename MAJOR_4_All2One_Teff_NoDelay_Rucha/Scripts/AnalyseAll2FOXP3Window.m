%% Analysing FOXP3 
% Intersection, applying new intesity filters, plot dynamics of
% selected upstream genes, plot biographs
% July 2017

%% FOXP3 PROBE SET 2 / 3    

GeneOfInterest = 'FOXP3';
IndexesOfProbessToUse = [2,3];
FileNameGeneOfInterestDynamics = [GeneOfInterest '_' CellType 'Cell12_Probes23.mat'];

for IndexGeneName = 1:size(IndexesOfProbessToUse,2)
    
    if IndexesOfProbessToUse(IndexGeneName) == 2
        CutoffGenesNumber = CutoffGenesNumberFOXP3_2_Window;
    elseif IndexesOfProbessToUse(IndexGeneName) == 3
        CutoffGenesNumber = CutoffGenesNumberFOXP3_3_Window;
    end

    %SelectedGeneName = GeneOfInterestProbesInCell12{IndexGeneName};
    VariableNameToLoad = [GeneOfInterest '_Probe' num2str(IndexesOfProbessToUse(IndexGeneName)) '_' CellType '12']; 
    SelectedGeneData = cell2mat(struct2cell(load(FileNameGeneOfInterestDynamics, VariableNameToLoad)));
    
%     % load file with resulting ranks for Cell 1, current probe
%     FileNameToLoad = [GeneOfInterest '_Results_Treg_AllGenes_Cell1ProbeSet_', ...
%         num2str(IndexesOfProbessToUse(IndexGeneName)),'.mat'];        
%     load(FileNameToLoad, 'ranks')
%     ranks2_Treg1 = ranks(1:CutoffGenesNumber,:);
 
    % load file with resulting ranks for Cell 2, current probe
    FileNameToLoad = [GeneOfInterest '_WINDOW_Results_' CellType '_AllGenes_Cell2ProbeSet_', ...
        num2str(IndexesOfProbessToUse(IndexGeneName)),'.mat'];        
    load(FileNameToLoad, 'ranks')
    ranks2_Treg2 = ranks(1:CutoffGenesNumber,:);

    % XXXXXX INTERSECT/UNION PROBES

    ProbeNumber = IndexesOfProbessToUse(IndexGeneName);

    % HERE WE DO THIS WITH ONLY CELL 2 !!!
    % Calling script to perform intersection and apply intensity filters
    [Filtered_CutIntersection_ProbeN,GeneNamesCell12, USELESS1, USELESS2,...
             ProbeIDsCell2, PublicIDsCell2] = IntersectSelectAnalysePlotProbesAll2oneWindow...
             (ranks2_Treg2, ranks2_Treg2,ThresholdAverageIntensityCell2,...
             ThresholdAverageIntensityCell2,CellType);
    
    % Calling script for LjungBox test, plot dynamics of filtered genes and Biographs
    load(['FILTERED_ExtractedCutData_' CellType 'ALLgen.mat'])
    if PlotBioGraph == 1   
        PlotDynamicsBiographsWindow(SelectedGeneData,Filtered_CutIntersection_ProbeN,...
         Tcell1DataTable,Tcell2DataTable,GeneNamesCell12, USELESS1, USELESS2,...
         ProbeIDsCell2, PublicIDsCell2,GeneOfInterest,ProbeNumber,PlotBioGraph,CellType)
    end
   
    if PlotDiGraph == 1
        ProbeNumber = [num2str(ProbeNumber) '-WINDOW'];
        plot_Digraph(Filtered_CutIntersection_ProbeN,GeneOfInterest,ProbeNumber,CellType);
    end

end



