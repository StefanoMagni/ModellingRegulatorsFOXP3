%% New code to first intersect then take union of Cell-1 and 2
% 27 Nov, 2017 (Rucha S. and Stefano M.)


clear; clc; close all;
load('FILTERED_ExtractedCutData_TregALLgen.mat')
load('FOXP3_Results_Treg_AllGenes_Cell1ProbeSet_2.mat', 'ranks')
RanksP2C1 = ranks;
load('FOXP3_Results_Treg_AllGenes_Cell1ProbeSet_3.mat', 'ranks')
RanksP3C1 = ranks;
load('FOXP3_Results_Treg_AllGenes_Cell2ProbeSet_2.mat', 'ranks')
RanksP2C2 = ranks;
load('FOXP3_Results_Treg_AllGenes_Cell2ProbeSet_3.mat', 'ranks')
RanksP3C2 = ranks;

% % clear ranks;
% % 
% % % Locations of genes in the main big data-file
% % indexP2C1 = cell2mat(RanksP2C1(:,1));
% % indexP3C1 = cell2mat(RanksP3C1(:,1));
% % indexP2C2 = cell2mat(RanksP2C2(:,1));
% % indexP3C2 = cell2mat(RanksP3C2(:,1));
% % 
% % %Extracting ProbeIDs & PublicIDs for the ranks above in cell 1 and 2
% % % Cell-1
% % ExtractProbeIDP2C1 = TableOfExtractedGeneNameAndIDs_cell1(indexP2C1,1:2);
% % ExtractProbeIDP3C1 = TableOfExtractedGeneNameAndIDs_cell1(indexP3C1,1:2);
% % % Cell-2
% % ExtractProbeIDP2C2 = TableOfExtractedGeneNameAndIDs_cell2(indexP2C2,1:2);
% % ExtractProbeIDP3C2 = TableOfExtractedGeneNameAndIDs_cell2(indexP3C2,1:2);
% % 
% % % Ranks with Probe ID and Public ID added in the last column no.5
% % Probe2Cell1All = [RanksP2C1(:,1:4) ExtractProbeIDP2C1 RanksP2C1(:,5)];
% % Probe3Cell1All = [RanksP3C1(:,1:4) ExtractProbeIDP3C1 RanksP3C1(:,5)];
% % Probe2Cell2All = [RanksP2C2(:,1:4) ExtractProbeIDP2C2 RanksP2C2(:,5)];
% % Probe3Cell2All = [RanksP3C2(:,1:4) ExtractProbeIDP3C2 RanksP3C2(:,5)];
% % 
% % % To create UNION of ProbeIDs from all 4 cases (P2C1, P3C1, P2C2, P3C2)
% % ProbeIDsP2C1 = Probe2Cell1All(:,5);
% % ProbeIDsP3C1 = Probe3Cell1All(:,5);
% % ProbeIDsP2C2 = Probe2Cell2All(:,5);
% % ProbeIDsP3C2 = Probe3Cell2All(:,5);
% % 
% % % Collect ProbeIDs that appear in at least one of: P2C1, P3C1, P2C2, P3C2
% % AllProbes = [ProbeIDsP2C1; ProbeIDsP3C1; ProbeIDsP2C2; ProbeIDsP3C2];
% % UniqueProbes = unique(AllProbes); 
% % Sets = {Probe2Cell1All, Probe3Cell1All, Probe2Cell2All, Probe3Cell2All};
% % 
% % 
% % %========================================
% % %% STEP 1:  copy info from each of 4 lists, put [] if no info
% % for jj = 1:size(Sets,2)
% %     MySet = Sets{jj};   
% %     
% %     for ii = 1:size(UniqueProbes,1)
% %         searchProbes = strcmp(UniqueProbes(ii,1),MySet(:,5));
% %         LocOfSearch = find(searchProbes == 1);
% %         
% %         if ~isempty(LocOfSearch)
% %             UniqueProbes{ii,1+jj} = MySet(LocOfSearch,:);
% %         else
% %             UniqueProbes{ii,1+jj} = []; %zeros(1,size(MySet,2)); 
% %         end
% %     end
% %     
% % end
% % 
% % 
% % %========================================
% % %% STEP 2:  discard lines which are not inside (P2c1 U P3c1) inters (P2c2 U P3c2)
% % %% STEP 3:  remember to apply activation/repression filter
% % UniqueProbesTag = [UniqueProbes, num2cell(ones(size(UniqueProbes,1),1))];
% % for kk = 1:size(UniqueProbes,1)
% %     CheckP2C1 = isempty(UniqueProbesTag{kk,2});
% %     CheckP3C1 = isempty(UniqueProbesTag{kk,3});
% %     CheckP2C2 = isempty(UniqueProbesTag{kk,4});
% %     CheckP3C2 = isempty(UniqueProbesTag{kk,5});   
% %     
% %     % P2, P3 both empty in cell indicated as '0'
% %     % TO CHECK THAT A PROBE APPEARS INTERSECTION OF BOTH CELLS
% %     if (CheckP2C1 == 1) && (CheckP3C1 == 1)
% %         UniqueProbesTag{kk,6} = 0;
% %     elseif (CheckP2C2 == 1) && (CheckP3C2 == 1)
% %         UniqueProbesTag{kk,6} = 0;
% %     end
% %     
% %     % CHECK IF ACT/REP IN ALL 4 CASES
% %     if UniqueProbesTag{kk,6} == 1
% %         CheckRepP2C1 = UniqueProbesTag{kk,2}{4};
% %         CheckRepP3C1 = UniqueProbesTag{kk,3}{4};
% %         CheckRepP2C2 = UniqueProbesTag{kk,4}{4};
% %         CheckRepP3C2 = UniqueProbesTag{kk,5}{4};
% %         if (CheckRepP2C1 >= 0)
% %             CheckRepP2C1 = 1; % positive
% %         elseif (CheckRepP2C1 < 0)
% %             CheckRepP2C1 = 0; % negative
% %         end
% %         if (CheckRepP3C1 >= 0)
% %             CheckRepP3C1 = 1; % positive
% %         elseif (CheckRepP3C1 < 0)
% %             CheckRepP3C1 = 0; % negative
% %         end
% %         if (CheckRepP2C2 >= 0)
% %             CheckRepP2C2 = 1; % positive
% %         elseif (CheckRepP2C2 < 0)
% %             CheckRepP2C2 = 0; % negative
% %         end
% %         if (CheckRepP3C2 >= 0)
% %             CheckRepP3C2 = 1; % positive
% %         elseif (CheckRepP3C2 < 0)
% %             CheckRepP3C2 = 0; % negative
% %         end 
% %         
% %         % For activation in all 4 sets Sum=4 / For repression in all 4 sets Sum=0
% %         Sum = CheckRepP2C1 + CheckRepP3C1 + CheckRepP2C2 + CheckRepP3C2;       
% %         if (Sum ~= 4) && (Sum ~= 0)
% %         	UniqueProbesTag{kk,6} = 0;
% %         end        
% %     end
% %     
% % end
% % 
% % Condition = (cell2mat(UniqueProbesTag(:,6)) == 0);  % Set a condition
% % UniqueProbesKeep = UniqueProbesTag(~Condition,:);   % Keep lines that meet the condition
% % UniqueProbesDiscard = UniqueProbesTag(Condition,:); % Remove lines that do not meet the condition
% %     
% % % Separating 4 cases and adding 'P2C1' etc. info in the last column
% % for ijkl = 1:size(UniqueProbesKeep,1)
% %     UniqueP2C1(ijkl,1:8) = [UniqueProbesKeep{ijkl,2}, 'P2-C1'];
% %     UniqueP3C1(ijkl,1:8) = [UniqueProbesKeep{ijkl,3}, 'P3-C1'];
% %     UniqueP2C2(ijkl,1:8) = [UniqueProbesKeep{ijkl,4}, 'P2-C2'];
% %     UniqueP3C2(ijkl,1:8) = [UniqueProbesKeep{ijkl,5}, 'P3-C2'];
% % end
% % 
% % 
% % filename = 'Array_P2P3_C1C2.mat';
% % save(filename)
% % save('Array_P2P3_C1C2.mat','UniqueP2C1','UniqueP3C1','UniqueP2C2','UniqueP3C2')
% % 
% % 
% % %========================================
% % %% STEP 4: order by fitness
% % 
% % CombineC1 = [UniqueP2C1; UniqueP3C1];  % Combine P2 and P3 in Cell-1
% % CombineC2 = [UniqueP2C2; UniqueP3C2];  % Combine P2 and P3 in Cell-2
% % % Combine P2 and P3 in Cell-1 and 2
% % InvertedCombineC1 = [UniqueP3C1; UniqueP2C1];  % Combine P2 and P3 in Cell-1
% % InvertedCombineC2 = [UniqueP3C2; UniqueP2C2];  % Combine P2 and P3 in Cell-2
% % CombineC1C2 = [CombineC1, CombineC2, InvertedCombineC1, InvertedCombineC2]; 
% % 
% % % RANKING CELL-2 (because selection is done for top candidates of Cell-2)
% % 
% % % Select 'cn': 1 for Cell-1 / 1+8 for Cell-2
% % cn = 1+8;
% % NewRank = [];
% % FileToLoad = CombineC1C2;
% % lineNumber = (1:size(FileToLoad,1))';
% % locationOfGenes = cell2mat(FileToLoad(lineNumber,cn));
% % fit_score = cell2mat(FileToLoad(lineNumber,cn+1));
% % names = FileToLoad(lineNumber,cn+2);
% % ActRep = cell2mat(FileToLoad(lineNumber,cn+3));
% % ProbeID1a = FileToLoad(lineNumber,cn+4);
% % RepPubID = FileToLoad(lineNumber,cn+5);
% % [ranks_Cell, idx] = ranking_list_FOXP3(fit_score,locationOfGenes,names,ActRep,ProbeID1a,RepPubID);
% % Ranks_Cell =  FileToLoad(idx,:);
% % 
% % Cell2_RankedForTable = [NewRank Ranks_Cell];   
% % Cell2_Ranked = Cell2_RankedForTable(:,1:16);
% % 
% % save('Cell2_RankedForTable.mat','Cell2_RankedForTable')
% % save('Cell2_Ranked.mat','Cell2_Ranked')


load('Cell2_RankedForTable.mat')

%========================================
%% STEP 5: CREATE TABLE FOR SUMMARISING INFO
% Tagging all Probes = 0, initially
Cell2_RankedForTable(:,33) = num2cell(zeros(size(Cell2_RankedForTable,1),1));

for jj = 1:size(Cell2_RankedForTable,1)
    if cell2mat(Cell2_RankedForTable(jj,33)) == 0   
        
        Cell2_RankedForTable{jj,33} = 1;
        
        GeneSearchInC1C2 = Cell2_RankedForTable(jj,5);
        SearchSameGene = strcmp(GeneSearchInC1C2,Cell2_RankedForTable(:,5));
        LocOfSAMEGene = find(SearchSameGene == 1);
        
        Cell2_RankedForTable{LocOfSAMEGene(2),33} = 99;      
    end
end

Condition = (cell2mat(Cell2_RankedForTable(:,33)) == 1);  % Set a condition
Cell2_ForTableKeep = Cell2_RankedForTable(Condition,:);   % Keep lines that meet the condition
Cell2_ForTableDiscard = Cell2_RankedForTable(~Condition,:);

for ijkl = 1:size(Cell2_ForTableKeep,1)
    ActRep = cell2mat(Cell2_ForTableKeep(ijkl,4));
    if ActRep >= 0
       Cell2_ForTableKeep{ijkl,33} = 'ACT';
    elseif ActRep < 0
       Cell2_ForTableKeep{ijkl,33} = 'REP';
    else
        disp('Error!')
    end
end

N = 250;
Cell2_TableToPrint = [Cell2_ForTableKeep(1:N,3), ...
                      Cell2_ForTableKeep(1:N,2), Cell2_ForTableKeep(1:N,5),...
                      Cell2_ForTableKeep(1:N,6), Cell2_ForTableKeep(1:N,8),...
                      Cell2_ForTableKeep(1:N,2+8), Cell2_ForTableKeep(1:N,5+8),...
                      Cell2_ForTableKeep(1:N,6+8), Cell2_ForTableKeep(1:N,8+8),... 
                      Cell2_ForTableKeep(1:N,2+16), Cell2_ForTableKeep(1:N,5+16),...
                      Cell2_ForTableKeep(1:N,6+16), Cell2_ForTableKeep(1:N,8+16),...   
                      Cell2_ForTableKeep(1:N,2+24), Cell2_ForTableKeep(1:N,5+24),...
                      Cell2_ForTableKeep(1:N,6+24), Cell2_ForTableKeep(1:N,8+24),...                     
                      Cell2_ForTableKeep(1:N,33)];

%========================================
%% STEP 6: SEARCH FOR TOP CANDIDATES AND PLOT

% Seach_All_Cases













