%% Analysing FOXP3 
% Intersection, applying new intesity filters, plot dynamics of
% selected upstream genes, plot biographs
% July 2017

%% FOXP3 PROBE SET 2 / 3    

GeneOfInterest = 'FOXP3';
IndexesOfProbessToUse = [2,3];
CutoffGenesNumber = CutoffGenesNumberFOXP3;
FileNameGeneOfInterestDynamics = [GeneOfInterest '_TregCell12_Probes23.mat'];

ranksUNION_Treg1 = [];
ranksUNION_Treg2 = [];
ProbeNumber = [];

for IndexGeneName = 1:size(IndexesOfProbessToUse,2)
        
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

    ranksUNION_Treg1 = [ranksUNION_Treg1; ranks2_Treg1];
    ranksUNION_Treg2 = [ranksUNION_Treg2; ranks2_Treg2];
    
    for i = 1:size(ranks2_Treg1)
        ProbeNumber = [ProbeNumber; IndexesOfProbessToUse(IndexGeneName), IndexesOfProbessToUse(IndexGeneName)];
    end
end

    % Calling script to perform intersection and apply intensity filters
    [Filtered_CutIntersection_ProbeN,GeneNamesCell12, ProbeIDsCell1, PublicIDsCell1,...
             ProbeIDsCell2, PublicIDsCell2, ProbeNumberCutOrdered] = ...
             IntersectSelectAnalysePlotProbesAll2oneUNION...
             (ranksUNION_Treg1, ranksUNION_Treg2, ThresholdAverageIntensityCell1,...
             ThresholdAverageIntensityCell2, ProbeNumber);
             
    % Calling script for LjungBox test, plot dynamics of filtered genes and Biographs
    load('FILTERED_ExtractedCutData_TregALLgen.mat')
    PlotDynamicsBiographsUNION(FileNameGeneOfInterestDynamics, Filtered_CutIntersection_ProbeN,...
         Tcell1DataTable,Tcell2DataTable,GeneNamesCell12, ProbeIDsCell1,...
         PublicIDsCell1,ProbeIDsCell2, PublicIDsCell2,GeneOfInterest,ProbeNumberCutOrdered,PlotBioGraph,SelectedGeneData)

    ProbeName = 's ALL';
    plot_Digraph(Filtered_CutIntersection_ProbeN,GeneOfInterest,ProbeName);




