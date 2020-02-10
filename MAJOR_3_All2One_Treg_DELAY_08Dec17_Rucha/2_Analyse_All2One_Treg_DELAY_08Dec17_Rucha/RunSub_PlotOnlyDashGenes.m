%% SEARCH AND PLOT GENES IN THE RANKED LIST OF Treg-1 and Treg-2
% 29 Nov, 2017 (Rucha S. and Stefano M.)

clc;
clear;
close all;

%---------------
STARTgene = 1;
ENDgene = 5;

SelectCandidates = ENDgene;

FromCell_1 = 0;
FromCell_2 = 1;

PlotFigure = 1;
%---------------
        
load('FILTERED_ExtractedCutData_TregALLgen.mat')
load('FOXP3_TregCell12_Probes23.mat')
load('CodingNonCodingList.mat');
load('Cell2_Ranked.mat')

CutOff = size(Cell2_Ranked,1);

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
    GenesToCheck = Cell1_Ranked(1:SelectCandidates,:);
end
if FromCell_2 == 1
    GenesToCheck = Cell2_Ranked(1:SelectCandidates,:);
end

% Creating a cutoff list for Cell 1,2
CutOffList = Cell2_Ranked(1:CutOff,:);

% Tagging 'SelectCandidates' = 0, initially
GenesToCheck(1:SelectCandidates,17) = num2cell(zeros(size(GenesToCheck,1),1));
 
%%===================================================
% Search and locate a each gene in cutoff Cell-1 list 

GeneSearchInC1C2 = '---';
SearchSameGene = strcmp(GeneSearchInC1C2,CutOffList(:,3));
LocOfSAMEGene = find(SearchSameGene == 1);

DashCutOff = [CutOffList(LocOfSAMEGene,:), num2cell(zeros(size(LocOfSAMEGene,1),1))];
        

for ii = ENDgene:-1:STARTgene
    
    if cell2mat(DashCutOff(ii,17)) == 0
        
        SearchDashProbe = strcmp(DashCutOff(ii,5),DashCutOff(:,5));
        LocOfDashProbe = find(SearchDashProbe == 1);
            
        % The genes is tagged as '1' at all places where it is located in the cutoff list
        for kl = 1:size(LocOfDashProbe,1)
            MyIndex_C1 = LocOfDashProbe(kl);
            DashCutOff{MyIndex_C1,17} = 1;
            DashCutOff{MyIndex_C1,18} = LocOfDashProbe;
        end
   

        %=======  CELL-1  =======
        % Collecting info and data-to-plot from the locations where the genes is found
        InfoOnPlot_C1 = DashCutOff(LocOfDashProbe,1:8);
        LocToPlot_C1 = cell2mat(DashCutOff(LocOfDashProbe,1));   % Location in the main data file
        DataToPlot_C1 = Tcell1DataTable(LocToPlot_C1,:);
        
        % Collecting info on ProbeID, PublicID, Fitness etc.
        Fitness_C1 = [repmat('Fit=',size(InfoOnPlot_C1,1),1),...
                      char(num2str(round(cell2mat(InfoOnPlot_C1(:,2)),1)))];
        ProbeIDs_C1 = [repmat(', ',size(InfoOnPlot_C1,1),1), char(strrep(InfoOnPlot_C1(:,5),'_','-'))];
        PublicIDs_C1 = [repmat(', ',size(InfoOnPlot_C1,1),1), char(strrep(InfoOnPlot_C1(:,6),'_','-'))];   
        CL_C1_round = round(LjungBoxTestQandCL_Tcell1(LocToPlot_C1,2),2);
        CL_C1 = [repmat(', CL=',size(InfoOnPlot_C1,1),1), char(num2str(CL_C1_round))];
        FOXP3probe_C1 = [repmat(', to ',size(InfoOnPlot_C1,1),1),...
                         num2str(cell2mat(InfoOnPlot_C1(:,8)))];
        for ijkl = 1:size(InfoOnPlot_C1,1)
            ActRep = cell2mat(InfoOnPlot_C1(ijkl,4));
            if ActRep >= 0
               InfoOnPlot_C1{ijkl,9} = 'ACT';
            elseif ActRep < 0
               InfoOnPlot_C1{ijkl,9} = 'REP';
            else
                disp('Error!')
            end
        end
        ActRepInfo_C1 = [repmat(', ',size(InfoOnPlot_C1,1),1), char(InfoOnPlot_C1(:,9))];
        
                     
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
        InfoOnPlot_C2 = DashCutOff(LocOfDashProbe,9:16);
        LocToPlot_C2 = cell2mat(DashCutOff(LocOfDashProbe,9));  % Location in the main data file
        DataToPlot_C2 = Tcell2DataTable(LocToPlot_C2,:);
                     
        % Collecting info on ProbeID, PublicID, Fitness etc.
        Fitness_C2 = [repmat('Fit=',size(InfoOnPlot_C2,1),1),...
                      char(num2str(round(cell2mat(InfoOnPlot_C2(:,2)),1)))];
        ProbeIDs_C2 = [repmat(', ',size(InfoOnPlot_C2,1),1), char(strrep(InfoOnPlot_C2(:,5),'_','-'))];
        PublicIDs_C2 = [repmat(', ',size(InfoOnPlot_C2,1),1), char(strrep(InfoOnPlot_C1(:,6),'_','-'))];   
        CL_C2_round = round(LjungBoxTestQandCL_Tcell2(LocToPlot_C2,2),2);
        CL_C2 = [repmat(', CL=',size(InfoOnPlot_C2,1),1), char(num2str(CL_C2_round))];
        FOXP3probe_C2 = [repmat(', to ',size(InfoOnPlot_C2,1),1),...
                         num2str(cell2mat(InfoOnPlot_C2(:,8)))];
        for pqr = 1:size(InfoOnPlot_C2,1)
            ActRep = cell2mat(InfoOnPlot_C2(pqr,4));
            if ActRep >= 0
               InfoOnPlot_C2{pqr,9} = 'ACT';
            elseif ActRep < 0
               InfoOnPlot_C2{pqr,9} = 'REP';
            else
                disp('Error!')
            end
        end
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

        
        
%% ===============      FIGURES      ==========================
        if PlotFigure == 1
            
            % CELL-1
            Legend_C1 = [Fitness_C1, FOXP3probe_C1, ProbeIDs_C1, PublicIDs_C1, CL_C1, ActRepInfo_C1];  

            figure;
%             h(1) = subplot(3,2,1);
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
%             h(3) = subplot(3,2,3);
            subplot(2,2,3)
            plot(time_points,GeneDataNORM_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
            hold on;
            plot(time_points,GeneDataNORM_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
            hold on;
            plot(time_points,y1NORM,'-*','Linewidth',1.8);
            xlim([0 360]); xlabel('time (min)'); ylabel('NORM. mRNA conc.');
            title(sprintf('%s', char(GeneSearchInC1C2), ' in ', num2str(cell2mat(InfoOnPlot_C1(1,8))),...
                          ' (NORMALISED)'),'FontSize',12);                                  
    %%===================================================    
        % CELL-2
            Legend_C2 = [Fitness_C2, FOXP3probe_C2, ProbeIDs_C2, PublicIDs_C2, CL_C2, ActRepInfo_C2];  

%             h(2) = subplot(3,2,2);
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
%             h(4) = subplot(3,2,4);
            subplot(2,2,4)
            plot(time_points,GeneDataNORM_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
            hold on;
            plot(time_points,GeneDataNORM_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
            hold on;
            plot(time_points,y2NORM,'-*','Linewidth',1.8);
            xlim([0 360]); xlabel('time (min)'); ylabel('NORM. mRNA conc.');
            title(sprintf('%s', char(GeneSearchInC1C2), ' in ', num2str(cell2mat(InfoOnPlot_C2(1,8))),...
                          ' (NORMALISED)'),'FontSize',12);
            
%             % Enlarge figure to full screen.
%             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0.5 0.8 0.8]);

    %%===================================================    
        % TABLE of PROBEID INFO

            CombinedInfo = [InfoOnPlot_C1, InfoOnPlot_C2];
            CombinedInfo(:,19) = num2cell(zeros(size(CombinedInfo,1),1));
            ProbeFitMatrix = [];
            ProbeNames = [];
            for rst = 1:size(CombinedInfo,1)
                if cell2mat(CombinedInfo(rst,19)) == 0
                    ProbeIDSearchFit = CombinedInfo(rst,5);
                    SearchProbeID = strcmp(ProbeIDSearchFit,CombinedInfo(:,5));
                    LocOfSAMEProbe = find(SearchProbeID == 1);
                    FitP2C1 = -1;
                    FitP3C1 = -1;
                    FitP2C2 = -1;
                    FitP3C2 = -1;
                    if (CombinedInfo{rst,8} == 'P2-C1')
                        FitP2C1 = CombinedInfo{rst,2};
                    elseif  (CombinedInfo{rst,8} == 'P3-C1')
                        FitP3C1 = CombinedInfo{rst,2};
                    end
                    if (CombinedInfo{rst,8+9} == 'P2-C2')
                        FitP2C2 = CombinedInfo{rst,2+9};
                    elseif  (CombinedInfo{rst,8+9} == 'P3-C2')
                        FitP3C2 = CombinedInfo{rst,2+9};
                    end
                    Loc2 = LocOfSAMEProbe(2);
                    if (CombinedInfo{Loc2,8} == 'P2-C1')
                        FitP2C1 = CombinedInfo{Loc2,2};
                    elseif  (CombinedInfo{Loc2,8} == 'P3-C1')
                        FitP3C1 = CombinedInfo{Loc2,2};
                    end
                    if (CombinedInfo{Loc2,8+9} == 'P2-C2')
                        FitP2C2 = CombinedInfo{Loc2,2+9};
                    elseif  (CombinedInfo{Loc2,8+9} == 'P3-C2')
                        FitP3C2 = CombinedInfo{Loc2,2+9};
                    end
                    MyVarName = {['Probe' ProbeIDSearchFit{1}]};
                    ProbeNames = [ProbeNames; MyVarName];     
                    ChangeName1 = strrep(ProbeNames(:,1),'_','');
                    ChangeName2 = strrep(ChangeName1,'-','');
                    ChangeName3 = strrep(ChangeName2,'/','');
                    ProbeN = [FitP2C1, FitP3C1, FitP2C2, FitP3C2]; 
                    ProbeFitMatrix = [ProbeFitMatrix; ProbeN];
                    CombinedInfo{rst,19} = 1;
                    CombinedInfo{Loc2,19} = 1;
                end
            end
% % 
% %             LastName = {'Cell1P2';'Cell1P3';'Cell2P2';'Cell2P3'};
% %             T = array2table(ProbeFitMatrix','VariableNames', ChangeName3);
% %             
% %             YourArray = table2array(T);
% %             YourNewTable = array2table(YourArray.');
% %             YourNewTable.Properties.RowNames = T.Properties.VariableNames;
% %           
% %             % PLOT TABLE in SUBPLOT
% % %             h(5) = subplot(3,2,5);
% %             uitable('Data',YourNewTable{:,:},'ColumnName',LastName,...
% %             'RowName',YourNewTable.Properties.RowNames,'Units', 'Normalized');% , 'Position',[0, 0, 1, 1]);% 
        end
    end
end  



%=========================================       
    