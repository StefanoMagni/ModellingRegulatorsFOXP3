%% Analysing CTLA4 
% Intersection, applying new intesity filters, plot dynamics of
% selected upstream genes, plot biographs
% July 2017

%% CTLA4 PROBE SET 1,2,3,5   

GeneOfInterest = 'CTLA4';
IndexesOfProbessToUse = [1,2,3,5];
FileNameGeneOfInterestDynamics = [GeneOfInterest '_TregCell12_Probes1235.mat'];

for IndexGeneName = 1:size(IndexesOfProbessToUse,2)
    
    if IndexesOfProbessToUse(IndexGeneName) == 1
        CutoffGenesNumber = CutoffGenesNumberCTLA4_1_Window;
    elseif IndexesOfProbessToUse(IndexGeneName) == 2
        CutoffGenesNumber = CutoffGenesNumberCTLA4_2_Window;
    elseif IndexesOfProbessToUse(IndexGeneName) == 3
        CutoffGenesNumber = CutoffGenesNumberCTLA4_3_Window;
    elseif IndexesOfProbessToUse(IndexGeneName) == 5 
        CutoffGenesNumber = CutoffGenesNumberCTLA4_5_Window;
    end

    %SelectedGeneName = GeneOfInterestProbesInCell12{IndexGeneName};
    VariableNameToLoad = [GeneOfInterest '_Probe' num2str(IndexesOfProbessToUse(IndexGeneName)) '_Treg12']; 
    SelectedGeneData = cell2mat(struct2cell(load(FileNameGeneOfInterestDynamics, VariableNameToLoad)));
    
%     % load file with resulting ranks for Cell 1, current probe
%     FileNameToLoad = [GeneOfInterest '_Results_Treg_AllGenes_Cell1ProbeSet_', ...
%         num2str(IndexesOfProbessToUse(IndexGeneName)),'.mat'];        
%     load(FileNameToLoad, 'ranks')
%     ranks2_Treg1 = ranks(1:CutoffGenesNumber,:);
 
    % load file with resulting ranks for Cell 2, current probe
    FileNameToLoad = [GeneOfInterest '_WINDOW_Results_Treg_AllGenes_Cell2ProbeSet_', ...
        num2str(IndexesOfProbessToUse(IndexGeneName)),'.mat'];        
    load(FileNameToLoad, 'ranks')
    ranks2_Treg2 = ranks(1:CutoffGenesNumber,:);

    ProbeNumber = IndexesOfProbessToUse(IndexGeneName);

    % HERE WE DO THIS WITH ONLY CELL 2 !!!
    % Calling script to perform intersection and apply intensity filters
    [Filtered_CutIntersection_ProbeN,GeneNamesCell12, USELESS1, USELESS2,...
             ProbeIDsCell2, PublicIDsCell2] = IntersectSelectAnalysePlotProbesAll2oneWindow...
             (ranks2_Treg2, ranks2_Treg2,ThresholdAverageIntensityCell2,...
             ThresholdAverageIntensityCell2);
    
    % Calling script for LjungBox test, plot dynamics of filtered genes and Biographs
    load('FILTERED_ExtractedCutData_TregALLgen.mat')
    PlotDynamicsBiographsWindow(SelectedGeneData,Filtered_CutIntersection_ProbeN,...
         Tcell1DataTable,Tcell2DataTable,GeneNamesCell12, USELESS1, USELESS2,...
         ProbeIDsCell2, PublicIDsCell2,GeneOfInterest,ProbeNumber,PlotBioGraph)
   
    ProbeNumber = [num2str(ProbeNumber) '-WINDOW'];
    plot_Digraph(Filtered_CutIntersection_ProbeN,GeneOfInterest,ProbeNumber);

end



