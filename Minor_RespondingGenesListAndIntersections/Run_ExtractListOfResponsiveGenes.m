%% Filtering

FilterData = 1;

ApplyFengFilters = 1;
ApplyAFFIMETRICflagsFilter = 1;
ApplyLjungBoxTestFilter = 1;
ApplyFilterToProbesWithMultipleGenes = 0;

% Feng Filters Parameters
ThrasholdFILTERING = 100;
NfoldFILTERING = 2; % EVENTUALLY NOT USED
ThrasholdMeanFILTERING = 50;
MinOverallDifferenceFILTERING = 100; % EVENTUALLY NOT USED
ThrasholdFengFILTERINGleftTail = 500; % EVENTUALLY NOT USED
% Ljung Box Test Parameters
LjungBoxTestCLTrashold = 0.99999; % 30.58 % Confidence Level for the Ljung Box Test, 1 means 100%
ApplyBonferroniCorrection = 1; % 0 = No, 1 = Yes.
Silent = 1; % 0 = No, 1 = Yes.

if FilterData
%     GeneDataFiltering('ExtractedCutData_Teff_KNOWN.mat', ApplyFengFilters, ...
%         ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
%         NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ...
%         ApplyFilterToProbesWithMultipleGenes, 1, 'Teff_KNOWN', ...
%         ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
%     GeneDataFiltering('ExtractedCutData_Treg_KNOWN.mat', ApplyFengFilters, ...
%         ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
%         NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ...
%         ApplyFilterToProbesWithMultipleGenes, 1, 'Treg_KNOWN', ...
%         ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
%     GeneDataFiltering('ExtractedCutData_Teff_MIGHT.mat', ApplyFengFilters, ...
%         ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
%         NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
%         ApplyFilterToProbesWithMultipleGenes, 1, 'Teff_MIGHT', ...
%         ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
%     GeneDataFiltering('ExtractedCutData_Treg_MIGHT.mat', ApplyFengFilters, ...
%         ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
%         NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
%         ApplyFilterToProbesWithMultipleGenes, 1, 'Treg_MIGHT', ...
%         ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
    GeneDataFiltering('ExtractedCutData_TeffALLgen.mat', ApplyFengFilters, ...
        ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
        NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
        ApplyFilterToProbesWithMultipleGenes, 1, 'TeffALLgen', ...
        ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
    GeneDataFiltering('ExtractedCutData_TregALLgen.mat', ApplyFengFilters, ...
        ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
        NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
        ApplyFilterToProbesWithMultipleGenes, 1, 'TregALLgen', ...
        ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
end

clear

%% ExtractingLists

% Teff
FileNameToBeLoaded = ['FILTERED_ExtractedCutData_TeffALLgen.mat']
load(FileNameToBeLoaded, 'TableOfExtractedGeneNameAndIDs_cell1')
load(FileNameToBeLoaded, 'TableOfExtractedGeneNameAndIDs_cell2')
load(FileNameToBeLoaded, 'Tcell1DataTable')
load(FileNameToBeLoaded, 'Tcell2DataTable')

[Intersection_Teff_Cell1Cell2, IndexesTeff1, IndexesTeff2] = intersect(TableOfExtractedGeneNameAndIDs_cell1(:,4), ...
    TableOfExtractedGeneNameAndIDs_cell2(:,4), 'stable');

disp(size(TableOfExtractedGeneNameAndIDs_cell1,1));
disp(size(TableOfExtractedGeneNameAndIDs_cell2,1));
disp(size(Intersection_Teff_Cell1Cell2,1));

Teff1_IntersectionGenes = Tcell1DataTable(IndexesTeff1,:);
Teff2_IntersectionGenes = Tcell2DataTable(IndexesTeff2,:);

% Treg
FileNameToBeLoaded = ['FILTERED_ExtractedCutData_TregALLgen.mat']
load(FileNameToBeLoaded, 'TableOfExtractedGeneNameAndIDs_cell1')
load(FileNameToBeLoaded, 'TableOfExtractedGeneNameAndIDs_cell2')
load(FileNameToBeLoaded, 'Tcell1DataTable')
load(FileNameToBeLoaded, 'Tcell2DataTable')

[Intersection_Treg_Cell1Cell2, IndexesTreg1, IndexesTreg2] = intersect(TableOfExtractedGeneNameAndIDs_cell1(:,4), ...
    TableOfExtractedGeneNameAndIDs_cell2(:,4), 'stable');
Treg1_IntersectionGenes = Tcell1DataTable(IndexesTreg1,:);
Treg2_IntersectionGenes = Tcell2DataTable(IndexesTreg2,:);

disp(size(TableOfExtractedGeneNameAndIDs_cell1,1));
disp(size(TableOfExtractedGeneNameAndIDs_cell2,1));
disp(size(Intersection_Treg_Cell1Cell2,1));

% Treg and Teff
[IntersectionResponsive_Treg_Teff, IndexesTeff, IndexesTreg] = intersect(Intersection_Teff_Cell1Cell2, Intersection_Treg_Cell1Cell2, 'stable');
disp(size(IntersectionResponsive_Treg_Teff,1));

load('ResponsiveGenesNames.mat')

Intersection111 = intersect(Teff_MIGHTresp_GenesNames, IntersectionResponsive_Treg_Teff, 'stable');
Intersection222 = intersect(Treg_MIGHTresp_GenesNames, IntersectionResponsive_Treg_Teff, 'stable');
disp(size((Intersection111),1))
disp(size((Intersection222),1))

Teff1_ALLIntersectionGenes = Teff1_IntersectionGenes(IndexesTeff,:);
Teff2_ALLIntersectionGenes = Teff2_IntersectionGenes(IndexesTeff,:);
Treg1_ALLIntersectionGenes = Treg1_IntersectionGenes(IndexesTreg,:);
Treg2_ALLIntersectionGenes = Treg2_IntersectionGenes(IndexesTreg,:);


%% Plot

% GenesDataExtractor(1, IntersectionResponsive_Treg_Teff, 'Treg', 1)

Tcell1DataTable = Teff1_ALLIntersectionGenes;
Tcell2DataTable = Teff2_ALLIntersectionGenes;

TimePointsList = [0.]; % (minutes)
for i = 1:(size(Tcell1DataTable,2)-1)
    TimePointsList = [ TimePointsList, TimePointsList(i) + 20.];
end

figure;
% Make subplot 1
subplot(2,1,1);
for i = 1:size(Tcell1DataTable,1)
    plot(TimePointsList,Tcell1DataTable(i,:),'-o'); hold on;
end
title('Cell 1'); xlabel('time (min)'); ylabel('mRNA concentration (a.u.?)');

% Make subplot 2
subplot(2,1,2)
for i = 1:size(Tcell2DataTable,1)
    plot(TimePointsList,Tcell2DataTable(i,:),'-o'); hold on;
end
title('Cell 2'); xlabel('time (min)'); ylabel('mRNA concentration (a.u.?)');
hold on;

% Make legend
if size(IntersectionResponsive_Treg_Teff,1) < 80
    ListOfLegendNames = [];
    for j = 1:size(IntersectionResponsive_Treg_Teff,1)
        ListOfLegendNames = [ListOfLegendNames, IntersectionResponsive_Treg_Teff(j,1)];
    end
    legend(ListOfLegendNames);
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeastoutside'); 
end
hold off;

% Make general plot title
annotation('textbox', [0 0.9 1 0.1], ...
'String', ['Gene expression time-series for two cells.'], ...
'EdgeColor', 'none', ...
'HorizontalAlignment', 'center');

% Give a proper name to the figure and save it
FigName = sprintf(['GenesExpressionINTERSECTIONTeff']);
savefig(FigName);

%%%%%%%%%

Tcell1DataTable = Treg1_ALLIntersectionGenes;
Tcell2DataTable = Treg2_ALLIntersectionGenes;

TimePointsList = [0.]; % (minutes)
for i = 1:(size(Tcell1DataTable,2)-1)
    TimePointsList = [ TimePointsList, TimePointsList(i) + 20.];
end

figure;
% Make subplot 1
subplot(2,1,1);
for i = 1:size(Tcell1DataTable,1)
    plot(TimePointsList,Tcell1DataTable(i,:),'-o'); hold on;
end
title('Cell 1'); xlabel('time (min)'); ylabel('mRNA concentration (a.u.?)');

% Make subplot 2
subplot(2,1,2)
for i = 1:size(Tcell2DataTable,1)
    plot(TimePointsList,Tcell2DataTable(i,:),'-o'); hold on;
end
title('Cell 2'); xlabel('time (min)'); ylabel('mRNA concentration (a.u.?)');
hold on;

% Make legend
if size(IntersectionResponsive_Treg_Teff,1) < 80
    ListOfLegendNames = [];
    for j = 1:size(IntersectionResponsive_Treg_Teff,1)
        ListOfLegendNames = [ListOfLegendNames, IntersectionResponsive_Treg_Teff(j,1)];
    end
    legend(ListOfLegendNames);
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeastoutside'); 
end
hold off;

% Make general plot title
annotation('textbox', [0 0.9 1 0.1], ...
'String', ['Gene expression time-series for two cells.'], ...
'EdgeColor', 'none', ...
'HorizontalAlignment', 'center');

% Give a proper name to the figure and save it
FigName = sprintf(['GenesExpressionINTERSECTIONTreg']);
savefig(FigName);


