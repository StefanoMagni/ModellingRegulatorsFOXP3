% Main script to run for ATA on Treg/Teff data from He 2012 
% Created by Stefano Magni on 21 February 2017, stefano.magni@uni.lu

clear; close all; clc;

%% Select what you would like to do
% 100000000 Just He 2012 comparison 
% 001000011 Not Filtered Network Inference and plots - KNOWN and MIGHT
% 001001011 FILTERED with FENG FILTERS - Network Inference and plots - KNOWN and MIGHT

% If you wish to scroll gene 1 by 1, run the "ScrollingGenesCL.m" script

% Which steps do you want to do?
ExtractAFFIMETRICflagsFromFiles = 1; %%%%%
PreliminaryComparisonWithPaper = 1; %%%%%
ExtractGenesData_KNOWN_MIGHT = 1; 
ExtractGenesData_ALL = 1; 
ComputeAndPlotLjungBoxTestStatisticsQandCL = 1;
FilterData = 1;
PlotLjungHistogramsAfterFiltering = 1;
ApplyAll2All = 1;
ApplyAll2AllON_FILTERED_DATA = 1;
PlotBioGraph = 1;
PlotBioGraphON_FILTERED_DATA = 1;

% Data Filtering Options
ApplyFengFilters = 1;
ApplyLjungBoxTestFilter = 0;
ApplyAFFIMETRICflagsFilter = 1;
ApplyFilterToProbesWithMultipleGenes = 0; % EVENTUALLY NOT TO BE USED

% Feng Filters Parameters
ThrasholdFILTERING = 100;
NfoldFILTERING = 2; % EVENTUALLY NOT USED
ThrasholdMeanFILTERING = 50;
MinOverallDifferenceFILTERING = 100; % EVENTUALLY NOT USED
ThrasholdFengFILTERINGleftTail = 500; % EVENTUALLY NOT USED
% Ljung Box Test Parameters
LjungBoxTestCLTrashold = 0.9999; % 30.58 % Confidence Level for the Ljung Box Test, 1 means 100%
ApplyBonferroniCorrection = 0; % 0 = No, 1 = Yes.
Silent = 1; % 0 = No, 1 = Yes.

% All2all Network Inference Options
method = 'ATA';
order = 1; %Order of system for ATA
init_opt = 'n4sid'; %Choice between {'iv','svf','gpmf','n4sid','all'}. See doc "tfestoptions"

% Network visual representation using biograph
threshold = 80; % Treshold for All to All fitness score
    
%% Extract Affi Flags
if ExtractAFFIMETRICflagsFromFiles
    run('FIlter_AFFIMETRIC_FLAGS.m');
end

%% Run comparison With Paper
if PreliminaryComparisonWithPaper
    run('ExtractHe2012plotsData.m');
end

%% Initialization
fprintf('Initialization...\n');

%
if ExtractGenesData_KNOWN_MIGHT
    % Extract Treg/Teff data according to 4 Feng's lists
    GenesDataExtractor(0, {}, '', 1)
end

% Extract All Data, needed also for Computing Ljung Test Statistics for ALL genes
if ExtractGenesData_ALL == 1
    % Teff
    GenesDataExtractor(2, {}, 'Teff', 0)
    % Treg
    GenesDataExtractor(2, {}, 'Treg', 0)
end

%% Compute and plot Ljung Box Test Statistics Q and associated Confidence Level (C.L.)
if ComputeAndPlotLjungBoxTestStatisticsQandCL == 1
    format long;
    run('LjungBoxTest.m');
    format short;
end

%% Filtering
if FilterData
    GeneDataFiltering('ExtractedCutData_Teff_KNOWN.mat', ApplyFengFilters, ...
        ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
        NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ...
        ApplyFilterToProbesWithMultipleGenes, 1, 'Teff_KNOWN', ...
        ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
    GeneDataFiltering('ExtractedCutData_Treg_KNOWN.mat', ApplyFengFilters, ...
        ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
        NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ...
        ApplyFilterToProbesWithMultipleGenes, 1, 'Treg_KNOWN', ...
        ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
    GeneDataFiltering('ExtractedCutData_Teff_MIGHT.mat', ApplyFengFilters, ...
        ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
        NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
        ApplyFilterToProbesWithMultipleGenes, 1, 'Teff_MIGHT', ...
        ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
    GeneDataFiltering('ExtractedCutData_Treg_MIGHT.mat', ApplyFengFilters, ...
        ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
        NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
        ApplyFilterToProbesWithMultipleGenes, 1, 'Treg_MIGHT', ...
        ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
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

%% Plot Ljung Box Test Statistics histograms after filtering to verify and for comparison
%disp('HERE!')
if PlotLjungHistogramsAfterFiltering == 1
    %disp('HERE HERE HERE')
    run('PlotLjungBoxHistogramsAfterFiltering.m');
end

%% System Identification
if ApplyAll2All | ApplyAll2AllON_FILTERED_DATA
    for KnownMight = {'KNOWN','MIGHT'}
        for TcellType = {'Teff','Treg'}
            fprintf(['...Starting ', TcellType{1}, ' ', KnownMight{1}, '...\n']);

            % Load data
            if ApplyAll2AllON_FILTERED_DATA == 1
                FileNameToBeLoaded = ['FILTERED_ExtractedCutData_' , TcellType{1}, '_', KnownMight{1}, '.mat'];
            elseif ApplyAll2AllON_FILTERED_DATA == 0
                FileNameToBeLoaded = ['ExtractedCutData_' , TcellType{1}, '_', KnownMight{1}, '.mat'];
            else
                disp('Error in FileterDataValue');
            end
            
            fprintf('Loading data...\n');
            load(FileNameToBeLoaded)

            %% Network inference (from Laurent's script)
            fprintf('Network inference...\n');
            opt = tfestOptions('InitMethod',init_opt); %Option for system ID ATA initial estimate technique. 

            for CellIndex = 1:2
                %Run network inference technique (Do not consider self regulation)
                                
                if CellIndex == 1
                    data = Tcell1DataTable;
                    %data = sscanf(sprintf('%s*', Tcell1DataTable), '%f*');
                elseif CellIndex == 2
                    data = Tcell2DataTable;
                end
                [fitness, models, IOs, AICs,dcgains] = just_tfest(order, 1, data, opt); %Only mRNA
                %%[fitness, models, IOs, AICs] = just_tfest(order, Ts, dataLI); %All species (proteins, mRNA)
                if ApplyAll2AllON_FILTERED_DATA == 0
                    mRNA_names = TableOfExtractedGeneNameAndIDs(:,4);%{'LHY mRNA';'TOC1 mRNA';'Y mRNA';'PRR9 mRNA';'PRR7 mRNA';'NI mRNA';'GI mRNA';'ZTL protein'}; 
                elseif ApplyAll2AllON_FILTERED_DATA == 1
                    if CellIndex == 1
                        mRNA_names = TableOfExtractedGeneNameAndIDs_cell1(:,4);
                    elseif CellIndex == 2
                        mRNA_names = TableOfExtractedGeneNameAndIDs_cell2(:,4);
                    end
                end
                ranks = ranking_list(fitness,mRNA_names,AICs);
                FileName = sprintf(['Results_', TcellType{1}, '_', KnownMight{1}, '_Cell' num2str(CellIndex) '.mat']);
                save(FileName, 'fitness', 'models', 'IOs', 'AICs', 'ranks', 'mRNA_names', 'dcgains');
                
%                 DCgainsName = sprintf(['tempDCgains', TcellType{1}, KnownMight{1}, '_Cell' num2str(CellIndex)]);
%                 varname = genvarname(DCgainsName)
%                 varname = dcgains
            end                           
        end
    end
end

if PlotBioGraph || PlotBioGraphON_FILTERED_DATA
    %% Plotting Biograph of the reconstructed networks
    fprintf('Plot Biograph of Reconstructed Network...\n');

    for KnownMight = {'KNOWN','MIGHT'}
        for TcellType = {'Teff','Treg'}
            fprintf(['...Starting ', TcellType{1}, ' ', KnownMight{1}, '...\n']);

            % Load desired dataset
            fprintf('Loading data...\n');
%             FileNameToBeLoaded = ['ExtractedCutData_', TcellType{1}, '_', KnownMight{1}, '.mat'];
%             load(FileNameToBeLoaded)
            % Load data
            if PlotBioGraphON_FILTERED_DATA == 1
                FileNameToBeLoaded = ['FILTERED_ExtractedCutData_' , TcellType{1}, '_', KnownMight{1}, '.mat'];
            elseif PlotBioGraphON_FILTERED_DATA == 0
                FileNameToBeLoaded = ['ExtractedCutData_' , TcellType{1}, '_', KnownMight{1}, '.mat'];
            else
                disp('Error in FileterDataValue');
            end
            load(FileNameToBeLoaded);
            
            for CellIndex = 1:2
%                 DCgainsName = sprintf(['tempDCgains', TcellType{1}, KnownMight{1}, '_Cell' num2str(CellIndex)]);
%                 varname = genvarname(DCgainsName)
%                 dcgains = varname              % Load desired results of Sys Id 
                if CellIndex == 1                    
                    ResultsFileNameToBeLoaded = ['Results_', TcellType{1}, '_', KnownMight{1}, '_Cell1.mat'];
                elseif CellIndex == 2
                    ResultsFileNameToBeLoaded = ['Results_', TcellType{1}, '_', KnownMight{1}, '_Cell2.mat'];
                end
                load(ResultsFileNameToBeLoaded);
                
                % Chose table according to filtering / no filtering
                if PlotBioGraphON_FILTERED_DATA == 0
                    MyTable = TableOfExtractedGeneNameAndIDs;
                elseif PlotBioGraphON_FILTERED_DATA == 1
                    if CellIndex == 1
                        MyTable = TableOfExtractedGeneNameAndIDs_cell1;
                    elseif CellIndex == 2
                        MyTable = TableOfExtractedGeneNameAndIDs_cell2;
                    end
                end
                
                % Format Genes Names for Plotting
                mRNAs_names = {MyTable{1,4}};
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
                end

                % Make a Bio Graph representing the reconstracted gene regulatory network
                fprintf('Plotting graph of reconstructed network...\n');
                
                if PlotBioGraphON_FILTERED_DATA == 0
                    FigName = sprintf(['GeneRegulatoryNetwork_', TcellType{1}, '_Cell', num2str(CellIndex), '_', KnownMight{1}]);
                elseif PlotBioGraphON_FILTERED_DATA == 1
                    FigName = sprintf(['GeneRegulatoryNetwork_', TcellType{1}, '_Cell', num2str(CellIndex), '_', KnownMight{1}, '_FILTERED']);
                end
                
                plot_biograph(fitness,threshold,mRNAs_names,FigName,dcgains);
            end
        end
    end
end




