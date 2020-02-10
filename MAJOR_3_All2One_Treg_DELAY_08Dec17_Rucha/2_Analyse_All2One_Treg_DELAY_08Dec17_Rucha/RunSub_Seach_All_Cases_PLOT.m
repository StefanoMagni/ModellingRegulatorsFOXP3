%% SEARCH AND PLOT INFO IN ALL 6-DELAY CASES (TREG-1&2) FOR MANUAL SELECTION
% Plot genes as per the ranked list of Treg-2
% 11 Dec, 2017 (By Rucha Sawlekar)

clc;
clear;
close all;

%---------------
%% SELECT NO. OF GENES TO PLOT AT A TIME

% Select START and END gene number 
% for eg: to plot gene no. 5 to 10, select STARTgene = 1 & ENDgene = 5;
STARTgene = 1;
ENDgene = 5;

PlotFigure = 1;      % To plot the gene dynamics
FlipRepressor = 0;   % To flip the normalised plot of a repressor gene
DelayHistogram = 1;  % To include delay-histogram as subplot with gene dynamics
DisplayDelay = 1;    % Display the delay in figure for which the gene has max. fit score.


ExcludeDashGenes = 1;  % EXCLUDES '---' GENES if set to 1
Exclude07Dec = 0;      % EXCLUDES Genes from previously graded list on 07Dec2017

SelectCandidates = 7000;
%---------------
        
load('FILTERED_ExtractedCutData_TregALLgen.mat')
load('FOXP3_TregCell12_Probes23.mat')
load('CodingNonCodingList.mat');
load('Cell2_Ranked.mat')
% load('RESULTS_13Dec_MOD.mat')
% MAIN FILE contaning info on max fit for each gene with how much delay, in Cell-1&2 
load('RESULTS_14Dec.mat')
load('Fitness_13TimePointsC1.mat')

% Loading fit.score etc. info after applying delays 0 to 100 min
load('Delay0_Cell1_FINAL.mat');   load('Delay0_Cell2_FINAL.mat'); 
load('Delay20_Cell1_FINAL.mat');  load('Delay20_Cell2_FINAL.mat');  
load('Delay40_Cell1_FINAL.mat');  load('Delay40_Cell2_FINAL.mat'); 
load('Delay60_Cell1_FINAL.mat');  load('Delay60_Cell2_FINAL.mat'); 
load('Delay80_Cell1_FINAL.mat');  load('Delay80_Cell2_FINAL.mat'); 
load('Delay100_Cell1_FINAL.mat'); load('Delay100_Cell2_FINAL.mat'); 
 

% Excluding '---' genes from last 30 genes (out of 7030 genes) as they are anyway 
% with 0 fit.score in C-2. Keeping other genes.
DiscardedRESULTS = RESULTS_14Dec(SelectCandidates+1:7030,:); 
SearchDashDisc = strcmp('---',DiscardedRESULTS(:,3)); 
LocNoDashDisc = find(SearchDashDisc ~= 1);
% Now creating a list of genes EXCLUDING '---' genes from last 30 genes of 'RESULTS_14Dec'
Last30WODash = DiscardedRESULTS(LocNoDashDisc,:);  

SelectedRESULTS = [RESULTS_14Dec(1:SelectCandidates,:); Last30WODash];

CodingNonCodSHORTlist = [];
time_points = 0:20:360;


%% ===================================================
% FOXP3 Probe-2,3 data in Cell-1,2
GeneData_C1 = [FOXP3_Probe2_Treg12(1,:);FOXP3_Probe3_Treg12(1,:)]; % FOXP3 Probe 2,3 in Cell-1
GeneData_C2 = [FOXP3_Probe2_Treg12(2,:);FOXP3_Probe3_Treg12(2,:)]; % FOXP3 Probe 2,3 in Cell-2

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


%% RANK FITNESS SCORE IN CELL-2 (FIT.SCORE USED FOR RANKING IS AFTER APPLYING DELAY)
fit_score = cell2mat(SelectedRESULTS(:,35));  % Column 35 refers to fit.score in C-2 after delay
names = SelectedRESULTS(:,11);
ActRep = SelectedRESULTS(:,12);
[ranks, idx] = ranking_list_FOXP3(fit_score,names,ActRep);

RESULTS_Ranked = SelectedRESULTS(idx,:);      % Ranked as per highest fitness in cell-2, AFTER delay
RESULTS = RESULTS_Ranked;                
ORIGINAL = Cell2_Ranked;                      % Ranked as per highest fitness in cell-2, BEOFRE delay

% Tagging 'SelectCandidates' = 0, initially
RESULTS(:,17) = num2cell(zeros(size(RESULTS,1),1));
ORIGINAL(:,17) = num2cell(zeros(size(ORIGINAL,1),1));



%% ---------------------
% Excluding Genes from previously graded list on 07Dec2017 (round-1 of grading genes)
if Exclude07Dec == 1
    load('GradingGenes_07Dec.mat')
    for gs = 1:size(GradingGenes_07Dec,1)
        SearchGraded = strcmp(GradingGenes_07Dec(gs,3),RESULTS(:,3));       
        LocGraded = find(SearchGraded == 1);
        for lg = 1:size(LocGraded,1)
            RESULTS{LocGraded(lg),17} = 1;
        end     
    end
end

% EXCLUDING '---' GENES
if ExcludeDashGenes == 1
    SearchDash = strcmp('---',RESULTS(:,3)); 
    LocNoDash = find(SearchDash ~= 1);
    RESULTS = RESULTS(LocNoDash,:);  % List of genes WITHOUT '---' genes
end

GenesToCheck = RESULTS(:,:);         % Genes are to be checked from this list named as 'GenesToCheck'

% EXTRA INFO (optional): To find number of common genes in first 200 of 'RESULTS' and graded list of 186 genes
LocOfGradedCommon = find(cell2mat(RESULTS(1:200,17)) == 1);
LocOfNewGenes = find(cell2mat(RESULTS(1:200,17)) == 0);

% ListNewToPrint = RESULTS(1:200,:);
             


%% -----------------------------------------
% MAIN FOR-LOOP TO SEARCH AND PLOT EACH GENE 

for ii = ENDgene:-1:STARTgene  %  STARTgene:ENDgene   %  
%     disp(ii)
    
    if cell2mat(RESULTS(ii,17)) == 0       
        %%===================================================
        % Search and locate a each gene in cutoff list (NOTE: gene names sequence in cell 1&2 is exactly the same)
        
        GeneSearchInC1C2 = GenesToCheck(ii,11);
        SearchSameGene = strcmp(GeneSearchInC1C2,RESULTS(:,3));
        LocOfSAMEGene = find(SearchSameGene == 1); 
        
        % The genes is tagged as '1' at all places in the cutoff list once it is located
        for klm = 1:size(LocOfSAMEGene,1)
            MyIndex_C1 = LocOfSAMEGene(klm);
            RESULTS{MyIndex_C1,17} = 1;
            ORIGINAL{ii,17} = 1;
            RESULTS{MyIndex_C1,37} = LocOfSAMEGene;
            
            ActRep = cell2mat(ORIGINAL(ii,12));
            if ActRep >= 0
               ORIGINAL{ii,18} = 'ACT';
            elseif ActRep < 0
               ORIGINAL{ii,18} = 'REP';
            end
        end
        
        MaxFitForDelay_C1 = char(RESULTS(LocOfSAMEGene(1),31));
        MaxFitForDelay_C2 = char(RESULTS(LocOfSAMEGene(1),34));
        
        %=======  CELL-1  =======
        % Collecting info and data-to-plot from the locations where the genes is found
        InfoOnPlot_C1 = RESULTS(LocOfSAMEGene,[1:8 18 31:32]);
        LocToPlot_C1 = cell2mat(RESULTS(LocOfSAMEGene,1));   % Location in the main data file 'Tcell1DataTable'
        DataToPlot_C1 = Tcell1DataTable(LocToPlot_C1,:);
        
        % This is just to avoid error in case the info on delay is missing (rare case but possible!)
        KeepLines_C1 = find(~cellfun(@isempty,InfoOnPlot_C1(:,10)));
        
        %=======  CELL-2  =======
        % Collecting info and data-to-plot from the locations where the genes is found
        InfoOnPlot_C2 = RESULTS(LocOfSAMEGene,[9:16 18 34:35]);
        LocToPlot_C2 = cell2mat(RESULTS(LocOfSAMEGene,9));   % Location in the main data file 'Tcell2DataTable'
        DataToPlot_C2 = Tcell2DataTable(LocToPlot_C2,:);
        
        KeepLines_C2 = find(~cellfun(@isempty,InfoOnPlot_C2(:,10)));       
        KeepLines = intersect(KeepLines_C1,KeepLines_C2);
        
        InfoOnPlot_C1 = InfoOnPlot_C1(KeepLines,:);
        InfoOnPlot_C2 = InfoOnPlot_C2(KeepLines,:);
        
        %=======  CELL-1  =======
        % Collecting info on ProbeID, PublicID, Fitness etc.
        Fitness_C1 = [repmat('Fit=',size(InfoOnPlot_C1,1),1),...
                      char(num2str(round(cell2mat(InfoOnPlot_C1(:,2)),1))),...
                      repmat(', New Fit-',size(InfoOnPlot_C1,1),1), char(InfoOnPlot_C1(:,10)),...
                      repmat('=',size(InfoOnPlot_C1,1),1),...
                      char(num2str(round(cell2mat(InfoOnPlot_C1(:,11)),1)))];            
        ProbeIDs_C1 = [repmat(', ',size(InfoOnPlot_C1,1),1), char(strrep(InfoOnPlot_C1(:,5),'_','-'))];
        PublicIDs_C1 = [repmat(', ',size(InfoOnPlot_C1,1),1), char(strrep(InfoOnPlot_C1(:,6),'_','-'))];   
        CL_C1_round = round(LjungBoxTestQandCL_Tcell1(LocToPlot_C1,2),2);
        CL_C1 = [repmat(', CL=',size(InfoOnPlot_C1,1),1), char(num2str(CL_C1_round))];
        FOXP3probe_C1 = [repmat(', to ',size(InfoOnPlot_C1,1),1),...
                         num2str(cell2mat(InfoOnPlot_C1(:,8)))];
        ActRepInfo_C1 = [repmat(', ',size(InfoOnPlot_C1,1),1), char(InfoOnPlot_C1(:,9))];  
%         V_Cell12(ii,:) = CutOffList(ii,[1:16 19]);
%         V_Cell12( all(cellfun(@isempty,V_Cell12),2), : ) = [];
                     
        % NORMALISE the upstream genes data
        count = 0;
        clear y1CurrentNORM y1NORM;
        for jj = 1:size(InfoOnPlot_C1,1)
            count = count+1;
            y1CurrentNORM(:,count) = (DataToPlot_C1(jj,:) - min(DataToPlot_C1(jj,:))) /...
                     (max(DataToPlot_C1(jj,:)) - min(DataToPlot_C1(jj,:)));
            y1NORM = y1CurrentNORM';
        end                   
        
%===================================================       
        %=======  CELL-2  =======
        % Collecting info on ProbeID, PublicID, Fitness etc.
        Fitness_C2 = [repmat('Fit=',size(InfoOnPlot_C2,1),1),...
                      char(num2str(round(cell2mat(InfoOnPlot_C2(:,2)),1))),...
                      repmat(', New Fit-',size(InfoOnPlot_C2,1),1), char(InfoOnPlot_C2(:,10)),...
                      repmat('=',size(InfoOnPlot_C2,1),1),...
                      char(num2str(round(cell2mat(InfoOnPlot_C2(:,11)),1)))];
        ProbeIDs_C2 = [repmat(', ',size(InfoOnPlot_C2,1),1), char(strrep(InfoOnPlot_C2(:,5),'_','-'))];
        PublicIDs_C2 = [repmat(', ',size(InfoOnPlot_C2,1),1), char(strrep(InfoOnPlot_C2(:,6),'_','-'))];   
        CL_C2_round = round(LjungBoxTestQandCL_Tcell2(LocToPlot_C2,2),2);
        CL_C2 = [repmat(', CL=',size(InfoOnPlot_C2,1),1), char(num2str(CL_C2_round))];
        FOXP3probe_C2 = [repmat(', to ',size(InfoOnPlot_C2,1),1),...
                         num2str(cell2mat(InfoOnPlot_C2(:,8)))];
        ActRepInfo_C2 = [repmat(', ',size(InfoOnPlot_C2,1),1), char(InfoOnPlot_C2(:,9))];
        
        % NORMALISE the upstream genes data
        count = 0;
        clear y2CurrentNORM y2NORM;
        for jj = 1:size(InfoOnPlot_C2,1)
            count = count+1;
            y2CurrentNORM(:,count) = (DataToPlot_C2(jj,:) - min(DataToPlot_C2(jj,:))) /...
                     (max(DataToPlot_C2(jj,:)) - min(DataToPlot_C2(jj,:)));
            y2NORM = y2CurrentNORM';
        end                            
                                   
        
%%===================================================
        % Find info on 'coding/non-coding/other'
        FindCodingNonCoding = strcmp(GeneSearchInC1C2{:},...
            CodingNonCodingList(:,2));
        locationOfFindResults = find(FindCodingNonCoding);
        CodingNonCodFLAGtemp = CodingNonCodingList(locationOfFindResults,:);
        % If no info is found on coding/non-coding then just put '-' instead
        if isempty(locationOfFindResults)
            CodingNonCodFLAGtemp = '-';
        end
        CodingNonCodFLAG = unique(CodingNonCodFLAGtemp(:,1));
        CodingNonCodSHORTlist = [CodingNonCodSHORTlist; CodingNonCodFLAG];   
        
        
%% PLOT FIGURES
        if PlotFigure == 1
           Grade = '';
           RankFitness_DELAY_HISTG_11Dec
        end
                          
    end
end

             
             
             
             