%% SEARCH FOR GENES IN THE RANKED LIST OF Treg-1 and Treg-2

clc;
clear;
close all;

%---------------
N = 10;
SelectCandidates = 21;
FromCell_1 = 0;
FromCell_2 = 1;

PlotFigure = 1;
%---------------
        
load('FILTERED_ExtractedCutData_TregALLgen.mat')
load('FOXP3_TregCell12_Probes23.mat')
load('Cell1_List.mat')
load('Cell2_List.mat')
load('CodingNonCodingList.mat');
% load('Final_P2P3_rank_Treg1.mat')
% load('Final_P2P3_rank_Treg2.mat')

CutOff_C1 = size(Cell1_List,1);
CutOff_C2 = size(Cell2_List,1);

CodingNonCodSHORTlist = [];
time_points = 0:20:360;


%%===================================================
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
if FromCell_1 == 1
    GenesToCheck = Cell1_List(1:SelectCandidates,:);
end
if FromCell_2 == 1
    GenesToCheck = Cell2_List(1:SelectCandidates,:);
end

% Creating a cutoff list for Cell 1,2
CutOffList_C1 = Cell1_List(1:CutOff_C1,:);
CutOffList_C2 = Cell2_List(1:CutOff_C2,:);
% Tagging 'SelectCandidates' = 0, initially
GenesToCheck(1:SelectCandidates,9) = num2cell(zeros(size(GenesToCheck,1),1));



for ii = SelectCandidates:-1:20 % ttt = N:SelectCandidates
%     ii = (SelectCandidates+N) - (ttt-1);
    
    if cell2mat(GenesToCheck(ii,9)) == 0
        
        %%===================================================
        % CELL-1 : Search and locate a single gene in cutoff Cell-1 list 

        GeneSearchInC1C2 = GenesToCheck(ii,3);
        SearchSameGene_C1 = strcmp(GeneSearchInC1C2,CutOffList_C1(:,3));
        LocOfSAMEGene_C1 = find(SearchSameGene_C1 == 1);
        V_Cell1{ii} = LocOfSAMEGene_C1; 
        
        % The genes is tagged as '1' at all places where it is located in the cutoff list
        for ijk = 1:size(LocOfSAMEGene_C1,1)
            MyIndex_C1 = LocOfSAMEGene_C1(ijk);
            CutOffList_C1{MyIndex_C1,9} = 1;
            GenesToCheck{MyIndex_C1,9} = 1;
            CutOffList_C1{MyIndex_C1,10} = LocOfSAMEGene_C1;
        end

        % Collecting info and data-to-plot from the locations where the genes is found
        InfoOnPlot_C1 = CutOffList_C1(LocOfSAMEGene_C1,:);
        LocToPlot_C1 = cell2mat(CutOffList_C1(LocOfSAMEGene_C1,1));   % Location in the main data file
        DataToPlot_C1 = Tcell1DataTable(LocToPlot_C1,:);
        
        % Collecting info on ProbeID, PublicID, Fitness etc.
        Fitness_C1 = [repmat('Fit=',size(InfoOnPlot_C1,1),1),...
                      char(num2str(round(cell2mat(InfoOnPlot_C1(:,2)),1)))];
        ProbeIDs_C1 = [repmat(', ',size(InfoOnPlot_C1,1),1), char(strrep(InfoOnPlot_C1(:,5),'_','-'))];
        PublicIDs_C1 = [repmat(', ',size(InfoOnPlot_C1,1),1), char(InfoOnPlot_C1(:,6))];   
        CL_C1_round = round(LjungBoxTestQandCL_Tcell1(LocToPlot_C1,2),2);
        CL_C1 = [repmat(', CL=',size(InfoOnPlot_C1,1),1), char(num2str(CL_C1_round))];
        FOXP3probe_C1 = [repmat(' to ',size(InfoOnPlot_C1,1),1),...
                         num2str(cell2mat(InfoOnPlot_C1(:,7)))];

        % NORMALISE the upstream genes data
        count = 0;
        clear y1CurrentNORM y1NORM;
        for jj = 1:size(InfoOnPlot_C1,1)
            count = count+1;
            y1CurrentNORM(:,count) = (DataToPlot_C1(jj,:) - min(DataToPlot_C1(jj,:))) /...
                     (max(DataToPlot_C1(jj,:)) - min(DataToPlot_C1(jj,:)));
            y1NORM = y1CurrentNORM';
        end
       
%%===================================================
        % CELL-2 : Search and locate a single gene in cutoff Cell-2 list 

        GeneSearchInC1C2 = GenesToCheck(ii,3);
        SearchSameGene_C2 = strcmp(GeneSearchInC1C2,CutOffList_C2(:,3));
        LocOfSAMEGene_C2 = find(SearchSameGene_C2 == 1);
        V_Cell2{ii} = LocOfSAMEGene_C2; 
        
        % The genes is tagged as '1' at all places where it is located in the cutoff list
        for lmn = 1:size(LocOfSAMEGene_C2,1)
            MyIndex_C2 = LocOfSAMEGene_C2(lmn);
            CutOffList_C2{MyIndex_C2,9} = 1;
            GenesToCheck{MyIndex_C2,9} = 1;
            CutOffList_C2{MyIndex_C2,10} = LocOfSAMEGene_C2;
        end

        % Collecting info and data-to-plot from the locations where the genes is found
        InfoOnPlot_C2 = CutOffList_C2(LocOfSAMEGene_C2,:);
        LocToPlot_C2 = cell2mat(CutOffList_C2(LocOfSAMEGene_C2,1));   % Location in the main data file
        DataToPlot_C2 = Tcell2DataTable(LocToPlot_C2,:);
        
        % Collecting info on ProbeID, PublicID, Fitness etc.
        Fitness_C2 = [repmat('Fit=',size(InfoOnPlot_C2,1),1),...
                      char(num2str(round(cell2mat(InfoOnPlot_C2(:,2)),1)))];
        ProbeIDs_C2 = [repmat(', ',size(InfoOnPlot_C2,1),1), char(strrep(InfoOnPlot_C2(:,5),'_','-'))];
%         ProbeIDs_C1 = [repmat(', ProbeID=',size(InfoOnPlot_C1,1),1), char(strrep(InfoOnPlot_C1(:,5),'_','-'))];
        PublicIDs_C2 = [repmat(', ',size(InfoOnPlot_C2,1),1), char(InfoOnPlot_C2(:,6))];   
        CL_C2_round = round(LjungBoxTestQandCL_Tcell2(LocToPlot_C2,2),2);
        CL_C2 = [repmat(', CL=',size(InfoOnPlot_C2,1),1), char(num2str(CL_C2_round))];
        FOXP3probe_C2 = [repmat(' to ',size(InfoOnPlot_C2,1),1),...
                         num2str(cell2mat(InfoOnPlot_C2(:,7)))];

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



%% ===============      FIGURES      ==========================
        if PlotFigure == 1
            if (~isempty(DataToPlot_C1))
                % CELL-1
                Legend_C1 = [Fitness_C1, FOXP3probe_C1, ProbeIDs_C1, PublicIDs_C1, CL_C1, ];  

                figure;
                subplot(2,2,1)
                % UNMORMALISED
                plot(time_points,GeneData_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,GeneData_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,DataToPlot_C1,'-*','Linewidth',1.8);
                maxInProbes_C1 = max([GeneData_C1(:);DataToPlot_C1(:)]);   
                xlim([0 360]); ylim([0 (maxInProbes_C1+maxInProbes_C1*(1/2))])
                xlabel('time (min)'); ylabel('mRNA conc.'); 
                title(sprintf('%s', char(GeneSearchInC1C2), ' in ', num2str(cell2mat(InfoOnPlot_C1(1,8)))),...
                              'FontSize',12);
                MyListForLegend_C1 = {'FOXP3 Probe-2', 'FOXP3 Probe-3'};
                for iii = 1:size(Legend_C1,1)
                    MyListForLegend_C1{iii+2} = Legend_C1(iii,:);
                end
                legend(MyListForLegend_C1, 'FontWeight','bold');

                % NORMALISED
                subplot(2,2,3)
                plot(time_points,GeneDataNORM_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,GeneDataNORM_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,y1NORM,'-*','Linewidth',1.8);
                xlim([0 360]); xlabel('time (min)'); ylabel('NORM. mRNA conc.');
                title(sprintf('%s', char(GeneSearchInC1C2), ' in ', num2str(cell2mat(InfoOnPlot_C1(1,8))),...
                              ' (NORMALISED)'),'FontSize',12);                      
            else
                figure;
                subplot(2,2,1)
                % UNMORMALISED
                plot(time_points,GeneData_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,GeneData_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
                hold on;
                plot([0 360],[500, 500],'k-','Linewidth',1.8);
                xlim([0 360]); ylim([0 600]);
                xlabel('time (min)'); ylabel('NORM. mRNA conc.');
                legend('FOXP3 Probe-2', 'FOXP3 Probe-3', 'No match!', 'FontWeight','bold');
                title([char(GeneSearchInC1C2), ', No match in Cell-1!'], 'FontSize',12)
                % NORMALISED
                subplot(2,2,3)
                plot(time_points,GeneDataNORM_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,GeneDataNORM_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
                hold on;
                plot([0 360],[0.6, 0.6],'k-','Linewidth',1.8);
                xlim([0 360]); xlabel('time (min)'); ylabel('NORM. mRNA conc.');
                title([char(GeneSearchInC1C2), ', No match in Cell-1! - (NORMALISED)'], 'FontSize',12)
            end
        %%===================================================    
            % CELL-2
            if (~isempty(DataToPlot_C2))
                Legend_C2 = [Fitness_C2, FOXP3probe_C2, ProbeIDs_C2, PublicIDs_C2, CL_C2, ];  

        %         figure;
                subplot(2,2,2)
                % UNMORMALISED
                plot(time_points,GeneData_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,GeneData_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,DataToPlot_C2,'-*','Linewidth',1.8);
                maxInProbes_C2 = max([GeneData_C2(:);DataToPlot_C2(:)]);   
                xlim([0 360]); ylim([0 (maxInProbes_C2+maxInProbes_C2*(1/2))])
                xlabel('time (min)'); ylabel('mRNA conc.'); 
                title(sprintf('%s', char(GeneSearchInC1C2), ' in ', num2str(cell2mat(InfoOnPlot_C2(1,8)))),...
                              'FontSize',12);
                MyListForLegend_C2 = {'FOXP3 Probe-2', 'FOXP3 Probe-3'};
                for jjj = 1:size(Legend_C2,1)
                    MyListForLegend_C2{jjj+2} = Legend_C2(jjj,:);
                end
                legend(MyListForLegend_C2, 'FontWeight','bold');

                % NORMALISED
                subplot(2,2,4)
                plot(time_points,GeneDataNORM_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,GeneDataNORM_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,y2NORM,'-*','Linewidth',1.8);
                xlim([0 360]); xlabel('time (min)'); ylabel('NORM. mRNA conc.');
                title(sprintf('%s', char(GeneSearchInC1C2), ' in ', num2str(cell2mat(InfoOnPlot_C2(1,8))),...
                              ' (NORMALISED)'),'FontSize',12);
            else
%                 figure;
                subplot(2,2,2)
                % UNMORMALISED
                plot(time_points,GeneData_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,GeneData_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
                hold on;
                plot([0 360],[500, 500],'k-','Linewidth',1.8);
                xlim([0 360]); ylim([0 600]);
                xlabel('time (min)'); ylabel('NORM.mRNA conc.');
                legend('FOXP3 Probe-2', 'FOXP3 Probe-3', 'No match!', 'FontWeight','bold');
                title([char(GeneSearchInC1C2), ', No match in Cell-2!'], 'FontSize',12)
                % NORMALISED
                subplot(2,2,4)
                plot(time_points,GeneDataNORM_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
                hold on;
                plot(time_points,GeneDataNORM_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
                hold on;
                plot([0 360],[0.6, 0.6],'k-','Linewidth',1.8);
                xlim([0 360]); xlabel('time (min)'); ylabel('NORM. mRNA conc.');
                title([char(GeneSearchInC1C2), ', No match in Cell-2! - (NORMALISED)'], 'FontSize',12)                      
            end
%             % Enlarge figure to full screen.
%             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0.5 0.8 0.8]);
        end
    end
end  



%=========================================       
    