%% Analysing IKZF4 
% Intersection, applying new intesity filters, plot dynamics of
% selected upstream genes, plot biographs
% July 2017

%% IKZF4 PROBE SET 2,3,4   

GeneOfInterest = 'IKZF4';
IndexesOfProbessToUse = [2,3,4];
FileNameGeneOfInterestDynamics = [GeneOfInterest '_TregCell12_Probes234.mat'];

for IndexGeneName = 1:size(IndexesOfProbessToUse,2)
    
    if IndexesOfProbessToUse(IndexGeneName) == 2
        CutoffGenesNumber = CutoffGenesNumberIKZF4_2;
    elseif IndexesOfProbessToUse(IndexGeneName) == 3
        CutoffGenesNumber = CutoffGenesNumberIKZF4_3;
    elseif IndexesOfProbessToUse(IndexGeneName) == 4 
        CutoffGenesNumber = CutoffGenesNumberIKZF4_4;
    end
    
    %SelectedGeneName = GeneOfInterestProbesInCell12{IndexGeneName};
    VariableNameToLoad = [GeneOfInterest '_Probe' num2str(IndexesOfProbessToUse(IndexGeneName)) '_Treg12']; 
    SelectedGeneData = cell2mat(struct2cell(load(FileNameGeneOfInterestDynamics, VariableNameToLoad)));
    
    % load file with resulting ranks for Cell 1, current probe
    FileNameToLoad = [GeneOfInterest '_Results_Treg_AllGenes_Cell1ProbeSet_', ...
        num2str(IndexesOfProbessToUse(IndexGeneName)),'.mat'];        
    load(FileNameToLoad, 'ranks')
    ranks2_Treg1 = ranks(1:CutoffGenesNumber,:);
 
    % load file with resulting ranks for Cell 2, current probe
    FileNameToLoad = [GeneOfInterest '_Results_Treg_AllGenes_Cell2ProbeSet_', ...
        num2str(IndexesOfProbessToUse(IndexGeneName)),'.mat'];        
    load(FileNameToLoad, 'ranks')
    ranks2_Treg2 = ranks(1:CutoffGenesNumber,:);

    % XXXXXX INTERSECT/UNION PROBES

    ProbeNumber = IndexesOfProbessToUse(IndexGeneName);

    % Calling script to perform intersection and apply intensity filters
    [Filtered_CutIntersection_ProbeN,GeneNamesCell12, ProbeIDsCell1, PublicIDsCell1,...
             ProbeIDsCell2, PublicIDsCell2] = IntersectSelectAnalysePlotProbesAll2one...
             (ranks2_Treg1, ranks2_Treg2,ThresholdAverageIntensityCell1,...
             ThresholdAverageIntensityCell2);
    
    % Calling script for LjungBox test, plot dynamics of filtered genes and Biographs
    load('FILTERED_ExtractedCutData_TregALLgen.mat')
    PlotDynamicsBiographs(SelectedGeneData,Filtered_CutIntersection_ProbeN,...
         Tcell1DataTable,Tcell2DataTable,GeneNamesCell12, ProbeIDsCell1,...
         PublicIDsCell1,ProbeIDsCell2, PublicIDsCell2,GeneOfInterest,ProbeNumber,PlotBioGraph)
   
    plot_Digraph(Filtered_CutIntersection_ProbeN,GeneOfInterest,ProbeNumber);
end



