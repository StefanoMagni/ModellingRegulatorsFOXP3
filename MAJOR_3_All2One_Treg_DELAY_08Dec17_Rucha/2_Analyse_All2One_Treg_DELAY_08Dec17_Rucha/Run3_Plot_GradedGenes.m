%% SEARCH AND PLOT GENES IN THE GRADED LIST OF Treg-1 and Treg-2
% 03 Dec, 2017 (By Rucha Sawlekar)

clc;
clear;
close all;

%==============

PlotFigure = 1;     % To plot figure
FlipRepressor = 0;  % To flip the normalised repressor gene plot for better decision
DelayHistogram = 1; % To include delay-histogram as subplot with gene dynamics
DelayAnalysis = 0;  % To plot no. of genes with highest fit.socre 
                    % in each delay case for 7000 & top 200 genes-list
DisplayDelay = 0;   % Display the delay in figure for which the gene has max. fit score.

% Grading done on 18Dec17 - FINAL (includes grades: 0, 100, 1, 1.5, 2, 8, 88, 0, 7777)
GradedGenes_18Dec = 1;  
% Grading done on 07Dec17 (includes grades: 0, 1, 1.5, 2, 2.5, 3, 3.5, 9, 99)
GradedGenes_07Dec = 0;

%==============

Category_100  = 1; % BEST in Both cells - STAR gene
Category_1  = 0; % 1 and 1L
Category_15 = 0; 
Category_2  = 0;
Category_25 = 0; % 2.5 and 2/3
Category_3  = 0;
Category_35 = 0; % 3Low and 3.5

Category_9  = 0; % SYNC-1
Category_99 = 0; % SYNC-2
Category_0  = 0; % REJECT

Category_8  = 0; % Trashed but Good in Cell-2
Category_88 = 0; % SYNC-2


%========================
if GradedGenes_18Dec == 1
    % Grading done on 18Dec17 - FINAL
    % (includes grades: 1, 1.5, 2, 2.5, 8, 88, 0, 7777)
    % Genes graded as '7777' means they are repeated
    load('GradedGenes_20Jan18.mat')
    GradingGenes = GradedGenes_20Jan18;
elseif GradedGenes_07Dec == 1
    % Grading done on 07Dec17
    % (includes grades: 1, 1.5, 2, 2.5, 3, 3.5, 9, 99, 0)
    load('GradingGenes_07Dec.mat')
    GradingGenes = GradingGenes_07Dec;
    DisplayDelay = 0;
end


%========================
% load('GradedGenes_18Dec_EDITED.mat')
% load('GradingGenes_07Dec.mat')
% load('GradingGenes_04Dec.mat')
% load('GradingGenes_05Dec.mat')
% load('GradingGenes_1.mat')
load('FILTERED_ExtractedCutData_TregALLgen.mat')
load('FOXP3_TregCell12_Probes23.mat')
load('CodingNonCodingList.mat');
load('Cell2_Ranked.mat')


CutOff = size(Cell2_Ranked,1);  % 7030 genes

CodingNonCodSHORTlist = [];
time_points = 0:20:360;

% Removing genes names showing '---' if any
SearchT = strcmp('---',GradingGenes(:,3));
LocOfT = find(SearchT == 1);
GradingGenes(LocOfT,:) = [];

%===================
% Selected and graded (14-17 Dec, 2017)
Loc1  = find(cell2mat(GradingGenes(:,18)) == 1);   % 1 and 1L
Loc15 = find(cell2mat(GradingGenes(:,18)) == 1.5);
Loc2  = find(cell2mat(GradingGenes(:,18)) == 2);
Loc25 = find(cell2mat(GradingGenes(:,18)) == 2.5);
Loc3  = find(cell2mat(GradingGenes(:,18)) == 3);
Loc35 = find(cell2mat(GradingGenes(:,18)) == 3.5);

% Rejected and Rejected-But-in-Sync (14-17 Dec, 2017)
Loc9  = find(cell2mat(GradingGenes(:,18)) == 9);   % SYNC - 1
Loc99 = find(cell2mat(GradingGenes(:,18)) == 99);  % SYNC - 2
Loc0  = find(cell2mat(GradingGenes(:,18)) == 0);   % REJECTED

% NEW GRADES FROM 18TH DEC 2017
Loc100  = find(cell2mat(GradingGenes(:,18)) == 100); % STAR GENE
Loc8  = find(cell2mat(GradingGenes(:,18)) == 8);     % TRASHED but good=1 ONLY in C-2
Loc88 = find(cell2mat(GradingGenes(:,18)) == 88);    % TRASHED but good=2 ONLY in C-2



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
%% ===================================================



% To display the 'Grade' in the title of figure which is given 
% to a gene (from top 200 list) after maunal selection based on dynamics
if Category_1 == 1
    GenesToCheck = GradingGenes(Loc1,:); Grade = 'Grade = 1';
elseif Category_15 == 1
    GenesToCheck = GradingGenes(Loc15,:); Grade = 'Grade = 1.5';
elseif Category_2 == 1
    GenesToCheck = GradingGenes(Loc2,:); Grade = 'Grade = 2';
elseif Category_25 == 1
    GenesToCheck = GradingGenes(Loc25,:); Grade = 'Grade = 2.5';
elseif Category_3 == 1
    GenesToCheck = GradingGenes(Loc3,:); Grade = 'Grade = 3';
elseif Category_35 == 1
    GenesToCheck = GradingGenes(Loc35,:); Grade = 'Grade = 3.5';
elseif Category_9 == 1
    GenesToCheck = GradingGenes(Loc9,:); Grade = 'Grade = SYNC-1'; 
elseif Category_99 == 1
    GenesToCheck = GradingGenes(Loc99,:); Grade = 'Grade = SYNC-2';  
elseif Category_0 == 1
    GenesToCheck = GradingGenes(Loc0,:); Grade = 'Grade = REJECT';  
    
elseif Category_100 == 1
    GenesToCheck = GradingGenes(Loc100,:); Grade = 'Grade = STAR GENE';     
elseif Category_8 == 1
    GenesToCheck = GradingGenes(Loc8,:); Grade = 'Grade = 1 in Only C2'; 
elseif Category_88 == 1
    GenesToCheck = GradingGenes(Loc88,:); Grade = 'Grade = 2 in Only C2';      
end

% Creating a cutoff list for Cell 1,2
CutOffList = Cell2_Ranked(1:CutOff,:);

% Tagging 'SelectCandidates' = 0, initially
GenesToCheck(:,18) = num2cell(zeros(size(GenesToCheck,1),1));

STARTgene = 1;
ENDgene = size(GenesToCheck,1);

for ii = ENDgene:-1:STARTgene  %  STARTgene:ENDgene  %   
    
    if cell2mat(GenesToCheck(ii,18)) == 0       
        %%===================================================
        % Search and locate a each gene in cutoff Cell-1 list 

        GeneSearchInC1C2 = GenesToCheck(ii,11);
        SearchSameGene = strcmp(GeneSearchInC1C2,CutOffList(:,3));
        % 'LocOfSAMEGene' gives row no.of gene in the list 'CutOffList' OR 'Cell2_Ranked', 
        % which is a list ranked according to Cell-2 ift. score & WITHOUT any delay
        LocOfSAMEGene = find(SearchSameGene == 1);  
        
        % The genes is tagged as '1' at all places where it is located in
        % the cutoff list -to avoide repetition
        for klm = 1:size(LocOfSAMEGene,1)
            MyIndex_C1 = LocOfSAMEGene(klm);
            CutOffList{MyIndex_C1,17} = 1;
            GenesToCheck{MyIndex_C1,18} = 1;
            CutOffList{MyIndex_C1,18} = LocOfSAMEGene;
            
            ActRep = cell2mat(CutOffList(LocOfSAMEGene,4));
            if ActRep(klm) >= 0
               CutOffList{MyIndex_C1,19} = 'ACT';
            elseif ActRep(klm) < 0
               CutOffList{MyIndex_C1,19} = 'REP';
            else
            end
        end

        if (Category_9 == 1) || (Category_99 == 1) || (Category_0 == 1)
            MaxFitForDelay_C1 = char(GenesToCheck{ii,2});
            MaxFitForDelay_C2 = char(GenesToCheck{ii,2+8});
        else
            MaxFitForDelay_C1 = char(GenesToCheck(ii,19));
            MaxFitForDelay_C2 = char(GenesToCheck(ii,21));
        end
        

        %=======  CELL-1  =======
        % Collecting info and data-to-plot from the locations where the genes is found
        InfoOnPlot_C1 = CutOffList(LocOfSAMEGene,[1:8 19]);
        LocToPlot_C1 = cell2mat(CutOffList(LocOfSAMEGene,1));   % Location in the main data file
        DataToPlot_C1 = Tcell1DataTable(LocToPlot_C1,:);
        
        % Collecting info on ProbeID, PublicID, Fitness etc.
        Fitness_C1 = [repmat('Fit=',size(InfoOnPlot_C1,1),1),...
                      char(num2str(round(cell2mat(InfoOnPlot_C1(:,2)),1)))];
        ProbeIDs_C1 = [repmat(', ProbeID=',size(InfoOnPlot_C1,1),1), char(strrep(InfoOnPlot_C1(:,5),'_','-'))];
        PublicIDs_C1 = [repmat(', PublicID=',size(InfoOnPlot_C1,1),1), char(strrep(InfoOnPlot_C1(:,6),'_','-'))];   
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
        % Collecting info and data-to-plot from the locations where the genes is found
        InfoOnPlot_C2 = CutOffList(LocOfSAMEGene,[9:16 19]);
        LocToPlot_C2 = cell2mat(CutOffList(LocOfSAMEGene,9));   % Location in the main data file
        DataToPlot_C2 = Tcell2DataTable(LocToPlot_C2,:);
                     
        % Collecting info on ProbeID, PublicID, Fitness etc.
        Fitness_C2 = [repmat('Fit=',size(InfoOnPlot_C2,1),1),...
                      char(num2str(round(cell2mat(InfoOnPlot_C2(:,2)),1)))];
        ProbeIDs_C2 = [repmat(', ProbeID=',size(InfoOnPlot_C2,1),1), char(strrep(InfoOnPlot_C2(:,5),'_','-'))];
        PublicIDs_C2 = [repmat(', PublicID=',size(InfoOnPlot_C2,1),1), char(strrep(InfoOnPlot_C1(:,6),'_','-'))];   
        CL_C2_round = round(LjungBoxTestQandCL_Tcell2(LocToPlot_C2,2),2);
        CL_C2 = [repmat(', CL=',size(InfoOnPlot_C2,1),1), char(num2str(CL_C2_round))];
        FOXP3probe_C2 = [repmat(', to ',size(InfoOnPlot_C2,1),1),...
                         num2str(cell2mat(InfoOnPlot_C2(:,8)))];
        ActRepInfo_C2 = [repmat(', ',size(InfoOnPlot_C1,1),1), char(InfoOnPlot_C2(:,9))];
                     
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



%% ===============      PLOT FIGURES      ==========================
        if PlotFigure == 1
           Grade = '';
           RankFitness_DELAY_HISTG_11Dec
           
       end
    end
end  


%% ===============      Plot Analysis     ================
if (PlotFigure == 1) && (DelayAnalysis == 1)
    
    load('RESULTS.mat')

%%  To plot no. of genes with highest fit.socre in each delay case for 7000 & top 200 gene-list
    delaystr = {'D-0', 'D-20', 'D-40', 'D-60', 'D-80', 'D-100'};
    delaystrSave = {'D_0', 'D_20', 'D_40', 'D_60', 'D_80', 'D_100'};  % Use '_' in the name to save mat file

    for DI = 1:size(delaystr,2)
        DelayCheck = delaystr{DI};

        for CellIndex = 1:2  % Cells => 2       
            if CellIndex == 1  
                col = 31;
                check = strcmp(DelayCheck,RESULTS(:,col));
                show = find(check == 1);               % Gives row no. of genes with the each delay
                show200 = find(check(1:200,1) == 1);   % Gives row no. of genes with the each delay within only top 200
                NoOfGenes_inTotal = size(show,1);
                NoOfGenes_inTop200 = size(show200,1);
                Total_inC1.(delaystrSave{DI}) = NoOfGenes_inTotal;   %for 7030 genes in C-1 (CELL format within STRUCTURE)
                Top200_inC1.(delaystrSave{DI}) = NoOfGenes_inTop200; %for 200 genes in C-1 (CELL format within STRUCTURE)
                Total_C1 = cell2mat(struct2cell(Total_inC1));        % Coverting from Structure >> Cell >> numeric
                Top200_C1 = cell2mat(struct2cell(Top200_inC1));
            end
            if CellIndex == 2
                col = 34;
                check = strcmp(DelayCheck,RESULTS(:,col));
                show = find(check == 1);                % Gives row no. of genes with the each delay
                show200 = find(check(1:200,1) == 1);    % Gives row no. of genes with the each delay within only top 200
                NoOfGenes_inTotal = size(show,1);
                NoOfGenes_inTop200 = size(show200,1);
                Total_inC2.(delaystrSave{DI}) = NoOfGenes_inTotal;   %for 7030 genes in C-2 (CELL format within STRUCTURE)
                Top200_inC2.(delaystrSave{DI}) = NoOfGenes_inTop200; %for 200 genes in C-2 (CELL format within STRUCTURE)
                Total_C2 = cell2mat(struct2cell(Total_inC2));        % Coverting from Structure >> Cell >> numeric
                Top200_C2 = cell2mat(struct2cell(Top200_inC2)); 
            end    
        end

    end

        
 % PLOT DELAY ANALYSIS
    Delay = [0; 20; 40; 60; 80; 100];
    hT_C1 = [];
    hT_C2 = [];
    hT_C3 = [];
    hT_C4 = [];

    figure;
    subplot(2,2,1)
    NoOfGenesC1 = num2str(Total_C1(:));
    b1 = bar(Delay,Total_C1,'FaceColor',[0 .5 .5],'EdgeColor',[0 .8 .8],'LineWidth',1.5);  % 'EdgeColor',[0 .9 .9]
    b1.FaceAlpha = 0.4; 
    hold on; grid on;
    xlabel('Delay (min)'); ylabel('No. of genes')
    hT_C1=[hT_C1,text(b1(:).XData+b1(:).XOffset,b1(:).YData,NoOfGenesC1, ...
          'VerticalAlignment','bottom','horizontalalign','center','FontWeight','bold', 'FontSize',13)];
    title({'No. of genes with highest fitness score ',...
           'for each delay case, among {\color[rgb]{.2 .4 1}all 7000 genes} in Cell-1'},...
           'FontWeight','bold', 'FontSize',13);
    % %     title(' {\color[rgb]{.2 .4 1} No. of genes for each delay among 7030 genes in Cell-1}');
    

    subplot(2,2,2)
    NoOfGenes200C1 = num2str(Top200_C1(:));
    b2 = bar(Delay,Top200_C1,'r','EdgeColor',[1 .4 .4],'LineWidth',1.5);
    b2.FaceAlpha = 0.4; 
    hold on; grid on;
    xlabel('Delay (min)'); ylabel('No. of genes')
    hT_C2=[hT_C2,text(b2(:).XData+b2(:).XOffset,b2(:).YData,NoOfGenes200C1, ...
          'VerticalAlignment','bottom','horizontalalign','center','FontWeight','bold', 'FontSize',13)];
    title({'No. of genes with highest fitness score ',...
           'for each delay case among {\color[rgb]{1 .2 .4}TOP-200 genes} in Cell-1'},...
           'FontWeight','bold', 'FontSize',13);


    subplot(2,2,3)
    NoOfGenes = num2str(Total_C2(:));
    b3 = bar(Delay,Total_C2,'FaceColor',[0 .5 .5],'EdgeColor',[0 .8 .8],'LineWidth',1.5);
    b3.FaceAlpha = 0.4; 
    hold on; grid on;
    xlabel('Delay (min)'); ylabel('No. of genes')
    hT_C3=[hT_C3,text(b3(:).XData+b3(:).XOffset,b3(:).YData,NoOfGenes, ...
          'VerticalAlignment','bottom','horizontalalign','center','FontWeight','bold', 'FontSize',13)];
    title({'No. of genes with highest fitness score ',...
           'for each delay case, among {\color[rgb]{.2 .4 1}all 7000 genes} in Cell-2'},...
           'FontWeight','bold', 'FontSize',13);


    subplot(2,2,4)
    NoOfGenes200 = num2str(Top200_C2(:));
    b4 = bar(Delay,Top200_C2,'r','EdgeColor',[1 .4 .4],'LineWidth',1.5);
    b4.FaceAlpha = 0.4; 
    hold on; grid on;
    xlabel('Delay (min)'); ylabel('No. of genes')
    hT_C4=[hT_C4,text(b4(:).XData+b4(:).XOffset,b4(:).YData,NoOfGenes200, ...
          'VerticalAlignment','bottom','horizontalalign','center','FontWeight','bold', 'FontSize',13)];
    title({'No. of genes with highest fitness score ',...
           'for each delay case among {\color[rgb]{1 .2 .4}TOP-200 genes} in Cell-2'},...
           'FontWeight','bold', 'FontSize',13);
    
    % To Plot enlarged figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0.5 0.75 0.9]);
end

%   *****************          END            ******************





% %             if (~isempty(DataToPlot_C1))  % (ActRepInfo_C1(1,1)) == '' %
% %                 % CELL-1
% %                 Legend_C1 = [Fitness_C1, FOXP3probe_C1, ProbeIDs_C1, PublicIDs_C1, CL_C1, ActRepInfo_C1];  
% % 
% %                 figure;
% %                 subplot(2,2,1)
% %                 % UNMORMALISED
% %                 plot(time_points,GeneData_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,GeneData_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,DataToPlot_C1,'-*','Linewidth',1.8);
% %                 maxInProbes_C1 = max([GeneData_C1(:);DataToPlot_C1(:)]);   
% %                 xlim([0 360]); ylim([0 (maxInProbes_C1+maxInProbes_C1*(1/2))])
% %                 xlabel('time (min)'); ylabel('mRNA conc.'); 
% %                 title(sprintf('%s', char(GeneSearchInC1C2), ' in ', num2str(cell2mat(InfoOnPlot_C1(1,8))),...
% %                               ' (',char(CodingNonCodFLAG),'), ', Grade), 'FontSize',12);
% %                 MyListForLegend_C1 = {'FOXP3 Probe-2', 'FOXP3 Probe-3'};
% %                 for iii = 1:size(Legend_C1,1)
% %                     MyListForLegend_C1{iii+2} = Legend_C1(iii,:);
% %                 end
% %                 legend(MyListForLegend_C1, 'FontWeight','bold', 'FontSize',11);
% % 
% %                 % NORMALISED
% %                 subplot(2,2,3)
% %                 plot(time_points,GeneDataNORM_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,GeneDataNORM_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,y1NORM,'-*','Linewidth',1.8);
% %                 xlim([0 360]); xlabel('time (min)'); ylabel('NORM. mRNA conc.');
% %                 title(sprintf('%s', char(GeneSearchInC1C2), ' in ', num2str(cell2mat(InfoOnPlot_C1(1,8))),...
% %                               ' (NORMALISED)'),'FontSize',12);                      
% %             else
% %                 figure;
% %                 subplot(2,2,1)
% %                 % UNMORMALISED
% %                 plot(time_points,GeneData_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,GeneData_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot([0 360],[500, 500],'k-','Linewidth',1.8);
% %                 xlim([0 360]); ylim([0 600]);
% %                 xlabel('time (min)'); ylabel('NORM. mRNA conc.');
% %                 legend('FOXP3 Probe-2', 'FOXP3 Probe-3', 'No match!', 'FontWeight','bold');
% %                 title([char(GeneSearchInC1C2), ', No match in Cell-1!'], 'FontSize',12)
% %                 % NORMALISED
% %                 subplot(2,2,3)
% %                 plot(time_points,GeneDataNORM_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,GeneDataNORM_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot([0 360],[0.6, 0.6],'k-','Linewidth',1.8);
% %                 xlim([0 360]); xlabel('time (min)'); ylabel('NORM. mRNA conc.');
% %                 title([char(GeneSearchInC1C2), ', No match in Cell-1! - (NORMALISED)'], 'FontSize',12)
% %             end
% %         %%===================================================    
% %             % CELL-2
% %             if (~isempty(DataToPlot_C2))
% %                 Legend_C2 = [Fitness_C2, FOXP3probe_C2, ProbeIDs_C2, PublicIDs_C2, CL_C2, ActRepInfo_C2];  
% % 
% %                 subplot(2,2,2)
% %                 % UNMORMALISED
% %                 plot(time_points,GeneData_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,GeneData_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,DataToPlot_C2,'-*','Linewidth',1.8);
% %                 maxInProbes_C2 = max([GeneData_C2(:);DataToPlot_C2(:)]);   
% %                 xlim([0 360]); ylim([0 (maxInProbes_C2+maxInProbes_C2*(1/2))])
% %                 xlabel('time (min)'); ylabel('mRNA conc.'); 
% %                 title(sprintf('%s', char(GeneSearchInC1C2), ' in ', num2str(cell2mat(InfoOnPlot_C2(1,8))),...
% %                               ' (',char(CodingNonCodFLAG),'), ', Grade), 'FontSize',12);
% %                 MyListForLegend_C2 = {'FOXP3 Probe-2', 'FOXP3 Probe-3'};
% %                 for jjj = 1:size(Legend_C2,1)
% %                     MyListForLegend_C2{jjj+2} = Legend_C2(jjj,:);
% %                 end
% %                 legend(MyListForLegend_C2, 'FontWeight','bold', 'FontSize',11);
% % 
% %                 % NORMALISED
% %                 subplot(2,2,4)
% %                 plot(time_points,GeneDataNORM_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,GeneDataNORM_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,y2NORM,'-*','Linewidth',1.8);
% %                 xlim([0 360]); xlabel('time (min)'); ylabel('NORM. mRNA conc.');
% %                 title(sprintf('%s', char(GeneSearchInC1C2), ' in ', num2str(cell2mat(InfoOnPlot_C2(1,8))),...
% %                               ' (NORMALISED)'),'FontSize',12);
% %             else
% %                 subplot(2,2,2)
% %                 % UNMORMALISED
% %                 plot(time_points,GeneData_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,GeneData_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot([0 360],[500, 500],'k-','Linewidth',1.8);
% %                 xlim([0 360]); ylim([0 600]);
% %                 xlabel('time (min)'); ylabel('NORM.mRNA conc.');
% %                 legend('FOXP3 Probe-2', 'FOXP3 Probe-3', 'No match!', 'FontWeight','bold');
% %                 title([char(GeneSearchInC1C2), ', No match in Cell-2!'], 'FontSize',12)
% %                 % NORMALISED
% %                 subplot(2,2,4)
% %                 plot(time_points,GeneDataNORM_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot(time_points,GeneDataNORM_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
% %                 hold on;
% %                 plot([0 360],[0.6, 0.6],'k-','Linewidth',1.8);
% %                 xlim([0 360]); xlabel('time (min)'); ylabel('NORM. mRNA conc.');
% %                 title([char(GeneSearchInC1C2), ', No match in Cell-2! - (NORMALISED)'], 'FontSize',12)                      
% %             end
% % %             % To Plot enlarged figure
% % %             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0.5 0.8 0.8]);







