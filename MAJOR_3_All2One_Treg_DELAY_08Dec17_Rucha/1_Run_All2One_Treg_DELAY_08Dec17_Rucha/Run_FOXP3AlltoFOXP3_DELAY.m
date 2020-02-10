%% TO RUN ALL2ONE ON 14K, 7K OR 200 GENES
% Modified by Rucha S. on Dec, 2017 for the new delay analysis
% Referred from the script by Stefano on July, 2017

close all; clear; clc;

% Filtering
FilterData = 0;
NetworkInference = 1;
ApplyRanking = 0;
PlotBioGraph = 0;

% Select cell type
Treg = 1;
Teff = 0;

% Select size of upstrem group
Genes_14k = 0;
Genes_7k = 0;
Genes_200 = 1;

% Select filters
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

%===========
if Treg == 1
    CellType = 'Treg';
end
if Teff == 1
    CellType = 'Teff';
end    
%===========

%% FILTER DATA (includes Bonferroni)
if FilterData == 1          
    KNOWNfileToLoad =   ['ExtractedCutData_' CellType '_KNOWN.mat'];
    MIGHTfileToLoad =   ['ExtractedCutData_' CellType '_MIGHT.mat'];
    ALLGENEfileToLoad = ['ExtractedCutData_' CellType 'ALLgen.mat'];
    CellTypeKNOWN =   [CellType '_KNOWN'];
    CellTypeMIGHT =   [CellType '_MIGHT'];
    CellTypeALLGENE = [CellType 'ALLgen'];
                
    GeneDataFiltering(KNOWNfileToLoad, ApplyFengFilters, ...
    ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
    NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ...
    ApplyFilterToProbesWithMultipleGenes, 1, CellTypeKNOWN, ...
    ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);

    GeneDataFiltering(MIGHTfileToLoad, ApplyFengFilters, ...
    ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
    NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
    ApplyFilterToProbesWithMultipleGenes, 1, CellTypeMIGHT, ...
    ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);

    GeneDataFiltering(ALLGENEfileToLoad, ApplyFengFilters, ...
    ApplyLjungBoxTestFilter, ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFILTERING, ...
    NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ... 
    ApplyFilterToProbesWithMultipleGenes, 1, CellTypeALLGENE, ...
    ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, 15);
end


%% Extracting Lists
threshold = 20;
genes = {'FOXP3_0', 'FOXP3_20', 'FOXP3_40', 'FOXP3_60', 'FOXP3_80', 'FOXP3_100'};

for IndexGeneName = 1:size(genes,2)
    
    SelectedGeneName = genes{IndexGeneName};   % Selects 'FOXP3_0', 'FOXP3_20' and so on

    % Calling script to extract desired portions of datafiles (gene expression data)
    % according to provided lists of gene names, and to data structure.
    GenesDataExtractor(1, {'FOXP3'}, CellType, 1)
    
    IndexesOfFOXP3sToUseCell1 = [2,3];
    IndexesOfFOXP3sToUseCell2 = [2,3];   
        
    if IndexGeneName == 1   % Delay = 0 minutes
       delay = 0;
    elseif IndexGeneName == 2   % Delay = 20 minutes
       delay = 1;
    elseif IndexGeneName == 3   % Delay = 40 minutes
       delay = 2;
    elseif IndexGeneName == 4   % Delay = 60 minutes
       delay = 3;
    elseif IndexGeneName == 5   % Delay = 80 minutes
       delay = 4;
    elseif IndexGeneName == 6   % Delay = 100 minutes
       delay = 5;
    else
        disp('Error')
    end
    
    % Everything of FOXP3, data, probeID etc. (old but IMP file)
    LoadExtractedFOXP3 = ['ExtractedCutData_' CellType 'Custom.mat'];
    load(LoadExtractedFOXP3)

    Tcell1DataTableFOXP3 = Tcell1DataTable;  % FOXP3 data from Cell-1
    Tcell2DataTableFOXP3 = Tcell2DataTable;  % FOXP3 data from Cell-2
    % ProbeID, PublicID etc. info
    TableOfExtractedGeneNameAndIDsFOXP3 = TableOfExtractedGeneNameAndIDs;
    save([SelectedGeneName 'data' CellType 's.mat'],'Tcell1DataTableFOXP3', 'Tcell2DataTableFOXP3', ...
        'TableOfExtractedGeneNameAndIDsFOXP3')

    %---------------------------------------------------
    % Loading FOXP3 data and the upstrem genes data here
    load([SelectedGeneName 'data' CellType 's.mat'])
    
    if Genes_14k == 1
        % CELL-1
        load(['FILTERED_ExtractedCutData_' CellType 'ALLgen.mat'], 'TableOfExtractedGeneNameAndIDs_cell1')
        load(['FILTERED_ExtractedCutData_' CellType 'ALLgen.mat'], 'Tcell1DataTable')
        % CELL-2
        load(['FILTERED_ExtractedCutData_' CellType 'ALLgen.mat'], 'TableOfExtractedGeneNameAndIDs_cell2')
        load(['FILTERED_ExtractedCutData_' CellType 'ALLgen.mat'], 'Tcell2DataTable')
        SaveName = '14kGenes';
    end
    if Genes_7k == 1
        % 3515 genes shortlisted from the union (Treg) of P2-C1, P3-C1, P2-C2, P3-C2
        load('Cell_1_INFO.mat') 
        load('Cell_2_INFO.mat')
        SaveName = '7kGenes';
    end
    if Genes_200 == 1
        % NOTE: THIS IS USEFUL TO RUN ONLY FOR TREG-1 because, 
        % Top 200 genes shortlisted from the union (Treg) of P2-C1, P3-C1, P2-C2, P3-C2
        % Here, code will anyway run for Treg-2 as well but only for 200 genes
        load('Genes200_C1.mat')
        load('Cell_2_INFO.mat')
        SaveName = '200Genes';
    end

    

    %% Network inference (from Laurent's script)
    fprintf('Network inference...\n');
    %opt = tfestOptions('InitMethod',init_opt); %Option for system ID ATA initial estimate technique. 

    if NetworkInference == 1
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
                % CELL-1
                if CellIndex == 1
                    if Genes_14k == 1     % for 14,000 genes
                        data = Tcell1DataTable;
                        mRNA_names = TableOfExtractedGeneNameAndIDs_cell1(:,4);
                        MyTestFOXP3 = Tcell1DataTableFOXP3(MultipleFOXP3sIndex,:);
                    elseif Genes_7k == 1  % for 7,000 genes
                        data = DataOfUpstream_C1;
                        mRNA_names = TableOfCell1(:,4);
                        MyTestFOXP3 = Tcell1DataTableFOXP3(MultipleFOXP3sIndex,:);
                    elseif Genes_200 == 1  % for 200 genes
                        data = DataInC1;
                        mRNA_names = InfoInC1(:,4);
                        MyTestFOXP3 = Tcell1DataTableFOXP3(MultipleFOXP3sIndex,1+CutN:19);
                    end
                 % CELL-2
                 elseif CellIndex == 2
                    if Genes_14k == 1     % for 14,000 genes
                        CutN = 0;
                        data = Tcell2DataTable;
                        mRNA_names = TableOfExtractedGeneNameAndIDs_cell2(:,4);
                        MyTestFOXP3 = Tcell2DataTableFOXP3(MultipleFOXP3sIndex,:);
                    elseif Genes_7k == 1  % for 7,000 genes
                        CutN = 0;
                        data = DataOfUpstream_C2;
                        mRNA_names = TableOfCell2(:,4);
                        MyTestFOXP3 = Tcell2DataTableFOXP3(MultipleFOXP3sIndex,:);
                    elseif Genes_200 == 1  % for 200 genes
                        data = DataOfUpstream_C2(1:200,:); % because we do not cut the time points like in C-1
                        mRNA_names = TableOfCell2(1:200,4);
                        MyTestFOXP3 = Tcell2DataTableFOXP3(MultipleFOXP3sIndex,1+CutN:19);
                    end     
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

                %%% RUN ALL2ONE HERE %%%
                %[fitness, models, IOs, AICs, dcgains] = improved_tfest(order, 1, inputs, outputs)
%                 disp(['Delay is = ' num2str(delay)])
                [fitness, models, IOs, AICs, dcgains] = all_to_one_delay_tfest_MOD11Dec(inputs,outputs,order,delay);
                % just_tfest(order, 1, data, opt); %Only mRNA

                MaxFit = max(fitness);              

                if ApplyRanking == 1
                    ranks = ranking_list_FOXP3(fitness,mRNA_names,AICs,dcgains);
                end       
                FileName = sprintf([SelectedGeneName '_' SaveName '_Results_' CellType '_DELAY_Cell'...
                           num2str(CellIndex) 'ProbeSet_' num2str(MultipleFOXP3sIndex) '.mat']);
                save(FileName, 'fitness', 'models', 'IOs', 'AICs', 'mRNA_names', 'dcgains','MaxFit');%, 'ranks');
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

                ResultsFileNameToBeLoaded = sprintf([SelectedGeneName '_' SaveName '_Results_' CellType '_DELAY_Cell'...
                           num2str(CellIndex) 'ProbeSet_' num2str(MultipleFOXP3sIndex) '.mat']);
                load(ResultsFileNameToBeLoaded); 

                if CellIndex == 1
                    MyTable = TableOfExtractedGeneNameAndIDs_cell1;
                elseif CellIndex == 2
                    MyTable = TableOfExtractedGeneNameAndIDs_cell2;
                end

                % Format Genes Names for Plotting
                mRNAs_names = {MyTable{1,4}};
                for i = 2:1:numel(MyTable(1:5,4))
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
                    SelectedGeneName '_' CellType '_Cell' num2str(CellIndex) 'ProbeSet_' ...
                    num2str(MultipleFOXP3sIndex)]);
                mRNA_name_SingleElement = SelectedGeneName;
                plot_biographMODforVectors(fitness,threshold,mRNAs_names,FigName,dcgains,mRNA_name_SingleElement); 

            end 
        end
    end
end


