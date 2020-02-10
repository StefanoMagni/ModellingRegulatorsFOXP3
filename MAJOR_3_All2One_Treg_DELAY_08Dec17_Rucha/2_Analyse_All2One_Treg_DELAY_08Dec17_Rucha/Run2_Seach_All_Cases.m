%% SEARCH FOLLOWING INFO IN ALL 6 CASES
% 11 Dec, 2017 (By Rucha Sawlekar)

% To search and collect information on: 
% [1] Fitness score of the gene for each of the 6 delay cases in Cell-1 and 2
% [2] Finding, storing and displaying the max fitness score among these 6
% delay cases in a single '.mat' file

% Weird genes locations [7007 7011 7013 7014 7018 7021 7023 7028 7029 7030]

clc;
clear;  
close all;
%---------------
STARTgene = 1;
ENDgene = 7030;

SelectCandidates = ENDgene;
%---------------
        
% Load files here
load('FILTERED_ExtractedCutData_TregALLgen.mat')
load('FOXP3_TregCell12_Probes23.mat')
load('CodingNonCodingList.mat');
load('Cell2_Ranked.mat')
% CutOff = size(Cell2_Ranked,1);

CodingNonCodSHORTlist = [];
time_points = 0:20:360;


%% ===================================================
% FOXP3 Probe-2,3 data in Cell-1,2
GeneData_C1 = [FOXP3_Probe2_Treg12(1,:);FOXP3_Probe3_Treg12(1,:)]; % Probe 2,3 in Cell-1
GeneData_C2 = [FOXP3_Probe2_Treg12(2,:);FOXP3_Probe3_Treg12(2,:)]; % Probe 2,3 in Cell-2

% Cell-1, Normlise FOXP3 probes
SelectedGeneNORMP2_C1 = (GeneData_C1(1,:) - min(GeneData_C1(1,:))) /...
        (max(GeneData_C1(1,:)) - min(GeneData_C1(1,:)));
SelectedGeneNORMP3_C1 = (GeneData_C1(2,:) - min(GeneData_C1(2,:))) /...
        (max(GeneData_C1(2,:)) - min(GeneData_C1(2,:)));
GeneDataNORM_C1 = [SelectedGeneNORMP2_C1;SelectedGeneNORMP3_C1];

% Cell-2, Normlise FOXP3 probes
SelectedGeneNORMP2_C2 = (GeneData_C2(1,:) - min(GeneData_C2(1,:))) /...
        (max(GeneData_C2(1,:)) - min(GeneData_C2(1,:)));
SelectedGeneNORMP3_C2 = (GeneData_C2(2,:) - min(GeneData_C2(2,:))) /...
        (max(GeneData_C2(2,:)) - min(GeneData_C2(2,:)));
GeneDataNORM_C2 = [SelectedGeneNORMP2_C2;SelectedGeneNORMP3_C2];   
%%===================================================


%% SELECT THE CELL FROM WHICH WE WANT TO FIND GENES

% 'Cell2_Ranked' is a list of 7030 genes (union of P2-C1, P3-C1, P2-C2, P3-C2)
% which are ranked in descending order of Cell-2 (P2&P3 both combined)fitness scores
GenesToCheck = Cell2_Ranked(1:SelectCandidates,:);


% Load delay cases results
load('Delay0_Cell1_FINAL.mat');   load('Delay0_Cell2_FINAL.mat'); 
load('Delay20_Cell1_FINAL.mat');  load('Delay20_Cell2_FINAL.mat');  
load('Delay40_Cell1_FINAL.mat');  load('Delay40_Cell2_FINAL.mat'); 
load('Delay60_Cell1_FINAL.mat');  load('Delay60_Cell2_FINAL.mat'); 
load('Delay80_Cell1_FINAL.mat');  load('Delay80_Cell2_FINAL.mat'); 
load('Delay100_Cell1_FINAL.mat'); load('Delay100_Cell2_FINAL.mat'); 

% Tagging 'SelectCandidates' = 0, initially
GenesToCheck(1:SelectCandidates,17) = num2cell(zeros(size(GenesToCheck,1),1));
             
for ii = STARTgene:ENDgene 
%     disp(ii)
    
    if cell2mat(GenesToCheck(ii,17)) == 0    % To avoid repetition 
        %%===================================================
        % Search and locate a each gene in cutoff Cell-1 list 

        GeneSearchInC1C2 = GenesToCheck(ii,11);
        SearchMain = strcmp(GeneSearchInC1C2,GenesToCheck(:,3));  % Searching for same 'Gene Name'
        LocInGenesToCheck = find(SearchMain == 1);
        
        % CELL-1
        SearchGene_0_C1 = strcmp(GeneSearchInC1C2,Delay0_Cell1_FINAL(:,2));
        LocOfGene_0_C1 = find(SearchGene_0_C1 == 1);
        
        SearchGene_20_C1 = strcmp(GeneSearchInC1C2,Delay20_Cell1_FINAL(:,2));
        LocOfGene_20_C1 = find(SearchGene_20_C1 == 1);
        
        SearchGene_40_C1 = strcmp(GeneSearchInC1C2,Delay40_Cell1_FINAL(:,2));
        LocOfGene_40_C1 = find(SearchGene_40_C1 == 1);
        
        SearchGene_60_C1 = strcmp(GeneSearchInC1C2,Delay60_Cell1_FINAL(:,2));
        LocOfGene_60_C1 = find(SearchGene_60_C1 == 1);
        
        SearchGene_80_C1 = strcmp(GeneSearchInC1C2,Delay80_Cell1_FINAL(:,2));
        LocOfGene_80_C1 = find(SearchGene_80_C1 == 1);
        
        SearchGene_100_C1 = strcmp(GeneSearchInC1C2,Delay100_Cell1_FINAL(:,2));
        LocOfGene_100_C1 = find(SearchGene_100_C1 == 1);            
        
        % The genes is tagged as '1' at all places where it is located in the cutoff list
        % CELL-1
        LocSet_C1 = {LocOfGene_0_C1, LocOfGene_20_C1, LocOfGene_40_C1,...
                  LocOfGene_60_C1, LocOfGene_80_C1, LocOfGene_100_C1};
        DelaySet_C1 = {Delay0_Cell1_FINAL, Delay20_Cell1_FINAL, Delay40_Cell1_FINAL,...
                  Delay60_Cell1_FINAL, Delay80_Cell1_FINAL, Delay100_Cell1_FINAL};
               
        count = 0;
        for iijj = 1:size(LocSet_C1,2) % 6 delay cases
            count = count+1;
            LocToCheck_C1 = LocSet_C1{iijj};
            SetToCheck_C1 = DelaySet_C1{iijj};           
            for abc = 1:size(LocToCheck_C1,1)  % no. of occurrences of selected gene   
                % size(LocSet_C2{1},1) 
                MyIndex_C1 = LocToCheck_C1(abc);
                GenesToCheck{LocInGenesToCheck(abc),17} = 1; 
                GenesToCheck{LocInGenesToCheck(abc),count+18} = SetToCheck_C1(LocToCheck_C1(abc),1);
%                 % To find ACT/REP
%                 ActRep = cell2mat(GenesToCheck(LocInGenesToCheck(abc),12));
%                 if ActRep >= 0
%                    GenesToCheck{LocInGenesToCheck(abc),18} = 'ACT';
%                 elseif ActRep < 0
%                    GenesToCheck{LocInGenesToCheck(abc),18} = 'REP';
%                 else
%                 end 
            end            
        end          

        
        % CELL-2
        SearchGene_0_C2 = strcmp(GeneSearchInC1C2,Delay0_Cell2_FINAL(:,2));
        LocOfGene_0_C2 = find(SearchGene_0_C2 == 1);
        
        SearchGene_20_C2 = strcmp(GeneSearchInC1C2,Delay20_Cell2_FINAL(:,2));
        LocOfGene_20_C2 = find(SearchGene_20_C2 == 1);
        
        SearchGene_40_C2 = strcmp(GeneSearchInC1C2,Delay40_Cell2_FINAL(:,2));
        LocOfGene_40_C2 = find(SearchGene_40_C2 == 1);
        
        SearchGene_60_C2 = strcmp(GeneSearchInC1C2,Delay60_Cell2_FINAL(:,2));
        LocOfGene_60_C2 = find(SearchGene_60_C2 == 1);
        
        SearchGene_80_C2 = strcmp(GeneSearchInC1C2,Delay80_Cell2_FINAL(:,2));
        LocOfGene_80_C2 = find(SearchGene_80_C2 == 1);
        
        SearchGene_100_C2 = strcmp(GeneSearchInC1C2,Delay100_Cell2_FINAL(:,2));
        LocOfGene_100_C2 = find(SearchGene_100_C2 == 1);
        
        % CELL-2
        LocSet_C2 = {LocOfGene_0_C2, LocOfGene_20_C2, LocOfGene_40_C2,...
                  LocOfGene_60_C2, LocOfGene_80_C2, LocOfGene_100_C2};
        DelaySet_C2 = {Delay0_Cell2_FINAL, Delay20_Cell2_FINAL, Delay40_Cell2_FINAL,...
                  Delay60_Cell2_FINAL, Delay80_Cell2_FINAL, Delay100_Cell2_FINAL};
               
        count = 0;
        for jjii = 1:size(LocSet_C2,2) % 6 delay cases
            count = count+1;
            LocToCheck_C2 = LocSet_C2{jjii};
            SetToCheck_C2 = DelaySet_C2{jjii};           
            for cba = 1:size(LocToCheck_C2,1)  % no. of occurrences of selected gene   
                % size(LocSet_C2{1},1) 
                MyIndex_C2 = LocToCheck_C2(cba); 
                GenesToCheck{LocInGenesToCheck(cba),count+24} = SetToCheck_C2(LocToCheck_C2(cba),1);
            end            
        end        
        
        FindMaxFit_C1 = GenesToCheck(LocInGenesToCheck,19:24);
        FindMaxFit_C2 = GenesToCheck(LocInGenesToCheck,25:30);
        % FIND MAX FITNESS AMONG ALL FITNESS SCORES
        for fm = 1:size(LocInGenesToCheck,1)
            % Creating array of all fitness scores
            FindMaxFitMat_C1 = cell2mat(cellfun(@(x) cell2mat(x),FindMaxFit_C1(fm,:),'un',0));
            FindMaxFitMat_C2 = cell2mat(cellfun(@(x) cell2mat(x),FindMaxFit_C2(fm,:),'un',0));
            % Find row & column indices of the maximum fitness score among above array
            [rowC1, colC1] = find(ismember(FindMaxFitMat_C1, max(FindMaxFitMat_C1(:))));
            [rowC2, colC2] = find(ismember(FindMaxFitMat_C2, max(FindMaxFitMat_C2(:))));
            % If all finess scores are '0' then to avoid error, row & column
            % index is selected as 1
            if all(FindMaxFitMat_C1 == 0)
               rowC1 = 1; colC1 = 1;
            elseif all(FindMaxFitMat_C2 == 0)
               rowC2 = 1; colC2 = 1;
            end
            % Add information for each gene, for eg., 'D-20' indicates,
            % max. fitness score for that gene is when delay is 20min
            % CELL-1
            if colC1 == 1
               GenesToCheck{LocInGenesToCheck(fm),31} = 'D-0'; 
               GenesToCheck{LocInGenesToCheck(fm),32} = max(FindMaxFitMat_C1(:));
               GenesToCheck{LocInGenesToCheck(fm),33} = {rowC1,colC1};
            elseif colC1 == 2
               GenesToCheck{LocInGenesToCheck(fm),31} = 'D-20';
               GenesToCheck{LocInGenesToCheck(fm),32} = max(FindMaxFitMat_C1(:));
               GenesToCheck{LocInGenesToCheck(fm),33} = {rowC1,colC1};
            elseif colC1 == 3
               GenesToCheck{LocInGenesToCheck(fm),31} = 'D-40';
               GenesToCheck{LocInGenesToCheck(fm),32} = max(FindMaxFitMat_C1(:));
               GenesToCheck{LocInGenesToCheck(fm),33} = {rowC1,colC1};
            elseif colC1 == 4
               GenesToCheck{LocInGenesToCheck(fm),31} = 'D-60';
               GenesToCheck{LocInGenesToCheck(fm),32} = max(FindMaxFitMat_C1(:));
               GenesToCheck{LocInGenesToCheck(fm),33} = {rowC1,colC1};
            elseif colC1 == 5
               GenesToCheck{LocInGenesToCheck(fm),31} = 'D-80';
               GenesToCheck{LocInGenesToCheck(fm),32} = max(FindMaxFitMat_C1(:));
               GenesToCheck{LocInGenesToCheck(fm),33} = {rowC1,colC1};
            elseif colC1 == 6
               GenesToCheck{LocInGenesToCheck(fm),31} = 'D-100';
               GenesToCheck{LocInGenesToCheck(fm),32} = max(FindMaxFitMat_C1(:));
               GenesToCheck{LocInGenesToCheck(fm),33} = {rowC1,colC1};
            else
            end
            
            % NOTE: 'ACT/REP' info is added using Cell-2 info, NOT using Cell-1
            % CELL-2
            if colC2 == 1
               GenesToCheck{LocInGenesToCheck(fm),34} = 'D-0';
               GenesToCheck{LocInGenesToCheck(fm),35} = max(FindMaxFitMat_C2(:));
               GenesToCheck{LocInGenesToCheck(fm),36} = {rowC2,colC2};
               ActRep = cell2mat(Delay0_Cell2_FINAL(LocOfGene_0_C2(1),3));
               if ActRep >= 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'ACT';
               elseif ActRep < 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'REP';
               else
               end
            elseif colC2 == 2
               GenesToCheck{LocInGenesToCheck(fm),34} = 'D-20';
               GenesToCheck{LocInGenesToCheck(fm),35} = max(FindMaxFitMat_C2(:));
               GenesToCheck{LocInGenesToCheck(fm),36} = {rowC2,colC2};
               ActRep = cell2mat(Delay20_Cell2_FINAL(LocOfGene_20_C2(1),3));
               if ActRep >= 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'ACT';
               elseif ActRep < 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'REP';
               else
               end
            elseif colC2 == 3
               GenesToCheck{LocInGenesToCheck(fm),34} = 'D-40';
               GenesToCheck{LocInGenesToCheck(fm),35} = max(FindMaxFitMat_C2(:));
               GenesToCheck{LocInGenesToCheck(fm),36} = {rowC2,colC2};
               ActRep = cell2mat(Delay40_Cell2_FINAL(LocOfGene_40_C2(1),3));
               if ActRep >= 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'ACT';
               elseif ActRep < 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'REP';
               else
               end
            elseif colC2 == 4
               GenesToCheck{LocInGenesToCheck(fm),34} = 'D-60';
               GenesToCheck{LocInGenesToCheck(fm),35} = max(FindMaxFitMat_C2(:));
               GenesToCheck{LocInGenesToCheck(fm),36} = {rowC2,colC2};
               ActRep = cell2mat(Delay60_Cell2_FINAL(LocOfGene_60_C2(1),3));
               if ActRep >= 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'ACT';
               elseif ActRep < 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'REP';
               else
               end
            elseif colC2 == 5
               GenesToCheck{LocInGenesToCheck(fm),34} = 'D-80';
               GenesToCheck{LocInGenesToCheck(fm),35} = max(FindMaxFitMat_C2(:));
               GenesToCheck{LocInGenesToCheck(fm),36} = {rowC2,colC2};
               ActRep = cell2mat(Delay80_Cell2_FINAL(LocOfGene_80_C2(1),3));
               if ActRep >= 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'ACT';
               elseif ActRep < 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'REP';
               else
               end
            elseif colC2 == 6
               GenesToCheck{LocInGenesToCheck(fm),34} = 'D-100';
               GenesToCheck{LocInGenesToCheck(fm),35} = max(FindMaxFitMat_C2(:));
               GenesToCheck{LocInGenesToCheck(fm),36} = {rowC2,colC2};
               ActRep = cell2mat(Delay100_Cell2_FINAL(LocOfGene_100_C2(1),3));
               if ActRep >= 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'ACT';
               elseif ActRep < 0
                   GenesToCheck{LocInGenesToCheck(fm),18} = 'REP';
               else
               end
            else
            end
        end
                       
    end
end








             