
close all; clear; clc;

% Filtering

FilterData = 1;
NetworkInference = 1;
ApplyRanking = 1;
PlotBioGraph = 0;

ApplyFengFilters = 1;
ApplyAFFIMETRICflagsFilter = 1;
ApplyLjungBoxTestFilter = 0;
ApplyFilterToProbesWithMultipleGenes = 0;

% Feng Filters Parameters
ThrasholdFILTERING = 100;
NfoldFILTERING = 2; % EVENTUALLY NOT USED
ThrasholdMeanFILTERING = 50;
MinOverallDifferenceFILTERING = 100; % EVENTUALLY NOT USED
ThrasholdFengFILTERINGleftTail = 500; % EVENTUALLY NOT USED
% Ljung Box Test Parameters
LjungBoxTestCLTrashold = 0.95; % 30.58 % Confidence Level for the Ljung Box Test, 1 means 100%
ApplyBonferroniCorrection = 0; % 0 = No, 1 = Yes.
Silent = 1; % 0 = No, 1 = Yes.

if FilterData
%     GeneDataFiltering('ExtractedCutData_Teff_KNOWN.mat', ApplyFengFilters, ...
%         ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
%         NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ...
%         ApplyFilterToProbesWithMultipleGenes, 1, 'Teff_KNOWN', ...
%         ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
    GeneDataFiltering('ExtractedCutData_Treg_KNOWN.mat', ApplyFengFilters, ...
        ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
        NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ...
        ApplyFilterToProbesWithMultipleGenes, 1, 'Treg_KNOWN', ...
        ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
%     GeneDataFiltering('ExtractedCutData_Teff_MIGHT.mat', ApplyFengFilters, ...
%         ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
%         NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
%         ApplyFilterToProbesWithMultipleGenes, 1, 'Teff_MIGHT', ...
%         ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
    GeneDataFiltering('ExtractedCutData_Treg_MIGHT.mat', ApplyFengFilters, ...
        ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
        NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
        ApplyFilterToProbesWithMultipleGenes, 1, 'Treg_MIGHT', ...
        ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
%     GeneDataFiltering('ExtractedCutData_TeffALLgen.mat', ApplyFengFilters, ...
%         ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
%         NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
%         ApplyFilterToProbesWithMultipleGenes, 1, 'TeffALLgen', ...
%         ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
    GeneDataFiltering('ExtractedCutData_TregALLgen.mat', ApplyFengFilters, ...
        ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
        NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
        ApplyFilterToProbesWithMultipleGenes, 1, 'TregALLgen', ...
        ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
end

%% ExtractingLists

threshold = 20;
genes = {'FOXP3', 'ICOS', 'CTLA4', 'IKZF4'};

for IndexGeneName = 1:size(genes,2)
    
    SelectedGeneName = genes{IndexGeneName};

    GenesDataExtractor(1, {SelectedGeneName}, 'Treg', 1)
    
    if IndexGeneName == 1   % FOXP3
        IndexesOfFOXP3sToUseCell1 = [2,3];
        IndexesOfFOXP3sToUseCell2 = [2,3];
    elseif IndexGeneName == 2  % ICOS
        IndexesOfFOXP3sToUseCell1 = [1];
        IndexesOfFOXP3sToUseCell2 = [1];
    elseif IndexGeneName == 3  % CTLA4
        IndexesOfFOXP3sToUseCell1 = [1,2,3,5];
        IndexesOfFOXP3sToUseCell2 = [1,2,3,5];
    elseif IndexGeneName == 4   % IKZF4
        IndexesOfFOXP3sToUseCell1 = [2,3,4];
        IndexesOfFOXP3sToUseCell2 = [2,3,4];    
    else
        disp('Error')
    end
    
    load('ExtractedCutData_TregCustom.mat')

    Tcell1DataTableFOXP3 = Tcell1DataTable
    Tcell2DataTableFOXP3 = Tcell2DataTable
    TableOfExtractedGeneNameAndIDsFOXP3 = TableOfExtractedGeneNameAndIDs
    save([SelectedGeneName 'dataTregs.mat'],'Tcell1DataTableFOXP3', 'Tcell2DataTableFOXP3', ...
        'TableOfExtractedGeneNameAndIDsFOXP3')

    load([SelectedGeneName 'dataTregs.mat'])
    % load('FILTERED_ExtractedCutData_TregALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell1')
    % load('FILTERED_ExtractedCutData_TregALLgen.mat', 'Tcell1DataTable')
    % load('FILTERED_ExtractedCutData_TregALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell2')
    % load('FILTERED_ExtractedCutData_TregALLgen.mat', 'Tcell2DataTable')
    load('FILTERED_ExtractedCutData_TregALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell1')
    load('FILTERED_ExtractedCutData_TregALLgen.mat', 'Tcell1DataTable')
    load('FILTERED_ExtractedCutData_TregALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell2')
    load('FILTERED_ExtractedCutData_TregALLgen.mat', 'Tcell2DataTable')

    %% Network inference (from Laurent's script)
    fprintf('Network inference...\n');
    %opt = tfestOptions('InitMethod',init_opt); %Option for system ID ATA initial estimate technique. 

    if NetworkInference
        for CellIndex = 1:2 %%%%%% FIX IT  %%%%%%
            % Run network inference technique (Do not consider self regulation)

            clear fitness models IOs AICs mRNA_names dcgains MaxFit ranks...
                  data dataNORM inputs;
            
            if CellIndex == 1
                IndexesOfFOXP3sToUse = IndexesOfFOXP3sToUseCell1;
            elseif CellIndex == 2
                IndexesOfFOXP3sToUse = IndexesOfFOXP3sToUseCell2;
            end

            for MultipleFOXP3sIndex = IndexesOfFOXP3sToUse
                if CellIndex == 1
                    data = Tcell1DataTable;%(1:100,:); % Tcell1DataTable; %%%%%% FIX IT  %%%%%%
                    mRNA_names = TableOfExtractedGeneNameAndIDs_cell1(:,4);
                    MyTestFOXP3 = Tcell1DataTableFOXP3(MultipleFOXP3sIndex,:);%+Tcell1DataTableFOXP3(3,:))/2; %%%%%% FIX IT  %%%%%%
                elseif CellIndex == 2
                    data = Tcell2DataTable;%(1:100,:); % Tcell2DataTable; %%%%%% FIX IT  %%%%%%
                    mRNA_names = TableOfExtractedGeneNameAndIDs_cell2(:,4);
                    MyTestFOXP3 = Tcell2DataTableFOXP3(MultipleFOXP3sIndex,:)%+Tcell2DataTableFOXP3(3,:))/2; %%%%%% FIX IT  %%%%%%
                end

                %%% Let's normalize %%%
                for i = 1:size(data,1)
                    dataNORM(i,:) = (data(i,:) - min(data(i,:))) / (max(data(i,:))- min(data(i,:)));
                end
                MyTestFOXP3NORM = (MyTestFOXP3 - min(MyTestFOXP3)) / (max(MyTestFOXP3) - min(MyTestFOXP3));
                %%% End Normalization %%%

                inputs = dataNORM;
                outputs = MyTestFOXP3NORM;
                order = 1;

                [fitness, models, IOs, AICs, dcgains] = improved_tfest(order, 1, inputs, outputs)
                % just_tfest(order, 1, data, opt); % Only mRNA

                MaxFit = max(fitness)

                % ranks = ranking_list(fitness,mRNA_names,AICs);
                if ApplyRanking == 1
                    ranks = ranking_list_FOXP3(fitness,mRNA_names,AICs,dcgains);
                end        
                FileName = sprintf([SelectedGeneName '_Results_Treg_AllGenes_Cell' num2str(CellIndex) 'ProbeSet_' num2str(MultipleFOXP3sIndex) '.mat']);
                save(FileName, 'fitness', 'models', 'IOs', 'AICs', 'mRNA_names', 'dcgains','MaxFit', 'ranks');
            end
        end   
    end


    if PlotBioGraph
        %% Plotting Biograph of the reconstructed networks
        fprintf('Plot Biograph of Reconstructed Network...\n');   
        for CellIndex = 1:2 

            if CellIndex == 1
                IndexesOfFOXP3sToUse = IndexesOfFOXP3sToUseCell1;
            elseif CellIndex == 2
                IndexesOfFOXP3sToUse = IndexesOfFOXP3sToUseCell2;
            end

            for MultipleFOXP3sIndex = IndexesOfFOXP3sToUse

                ResultsFileNameToBeLoaded = sprintf([SelectedGeneName '_Results_Treg_AllGenes_Cell' ...
                    num2str(CellIndex) 'ProbeSet_' num2str(MultipleFOXP3sIndex) '.mat']);
                load(ResultsFileNameToBeLoaded); 

                if CellIndex == 1
                    MyTable = TableOfExtractedGeneNameAndIDs_cell1;
                elseif CellIndex == 2
                    MyTable = TableOfExtractedGeneNameAndIDs_cell2;
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

                % Make a Bio Graph representing the reconstructed gene regulatory network
                fprintf('Plotting graph of reconstructed network...\n');

                FigName = sprintf([SelectedGeneName 'GeneRegulatoryNetwork_' ...
                    SelectedGeneName '_Treg_Cell' num2str(CellIndex) 'ProbeSet_' ...
                    num2str(MultipleFOXP3sIndex)]);
                mRNA_name_SingleElement = SelectedGeneName;
                plot_biographMODforVectors(fitness,threshold,mRNAs_names,FigName,dcgains,mRNA_name_SingleElement); 

            end 
        end
    end
end
% % % ranks = ranking_list_FOXP3(fitness,mRNA_names,AICs,dcgains);
% % % compare(IOs{10921},models{10921})
% % % figure; plot(IOs{10921})
% % 
% % % DataPointsVector = Tcell2DataTableFOXP3(3,:);
% % % 
% % % TestStatistics = ComputeLjungBoxTestStatistics(DataPointsVector, 15);
% % % %ListOfTestStatisticValues = [ListOfTestStatisticValues, TestStatistics]
% % % 
% % % ConfidenceLevels = chi2cdf(TestStatistics,15)
% % %             
% % % %%%%%%%%%% DEFINE FUNCTION TO COMPUTE LJUNG BOX TEST STATISTIC %%%%%%%%%%
% % % function [TestStatistics] = ComputeLjungBoxTestStatistics(DataPointsVector, NLags)
% % % Ndatapoins = numel(DataPointsVector);
% % % MyGeneAutocorrelation = autocorr(DataPointsVector); 
% % % TestStatistics = Ndatapoins * (Ndatapoins+2) * ...
% % %     sum(MyGeneAutocorrelation(2:NLags+1).^2 ./ ...
% % %     ((Ndatapoins-1):-1:(Ndatapoins-NLags)));
% % % end
