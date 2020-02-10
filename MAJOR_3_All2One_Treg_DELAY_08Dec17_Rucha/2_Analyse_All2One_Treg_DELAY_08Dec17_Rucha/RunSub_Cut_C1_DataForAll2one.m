% COllect data of 10 genes for trial

clc;
clear;

CutN = 6;

load('FILTERED_ExtractedCutData_TregALLgen.mat')
load('FOXP3_TregCell12_Probes23.mat')
load('Check_C1_All2One.mat')
load('RESULTS.mat')

% UPSTREAM dash genes



% UPSTREAM
LocationOf200Genes = cell2mat(RESULTS(1:200,1));
LocationOf700Genes = cell2mat(RESULTS(:,1));

% % 200 genes
% InfoInC1 = TableOfExtractedGeneNameAndIDs_cell1(LocationOf200Genes,:);
% DataInC1 = Tcell1DataTable(LocationOf200Genes,1+CutN:19); % 13 data points (7 to 19)

% 7000 genes
InfoInC1 = TableOfExtractedGeneNameAndIDs_cell1(LocationOf700Genes,:);
DataInC1 = Tcell1DataTable(LocationOf700Genes,1+CutN:19); % 13 data points (7 to 19)


% % FOXP3
% GeneData_C1 = [FOXP3_Probe2_Treg12(1,1+CutN:19);FOXP3_Probe3_Treg12(1,1+CutN:19)]; % Probe 2,3 in Cell-1

% save('Genes200_C1.mat','InfoInC1','DataInC1','CutN')
save('Genes7000_C1.mat','InfoInC1','DataInC1','CutN')






    
    
    