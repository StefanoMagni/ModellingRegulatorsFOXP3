%% %%  Plot figures and histogram for DELAYS
% 11 Dec, 2017 (By Rucha Sawlekar)

Delay = [0 20 40 60 80 100];

if DelayHistogram == 1
% Creating a cutoff list for Cell 1,2
    load('Delay0_Cell1_FINAL.mat');   load('Delay0_Cell2_FINAL.mat'); 
    load('Delay20_Cell1_FINAL.mat');  load('Delay20_Cell2_FINAL.mat');  
    load('Delay40_Cell1_FINAL.mat');  load('Delay40_Cell2_FINAL.mat'); 
    load('Delay60_Cell1_FINAL.mat');  load('Delay60_Cell2_FINAL.mat'); 
    load('Delay80_Cell1_FINAL.mat');  load('Delay80_Cell2_FINAL.mat'); 
    load('Delay100_Cell1_FINAL.mat'); load('Delay100_Cell2_FINAL.mat'); 

    %%===================================================
    % Search and locate a each gene in cutoff Cell-1 list     
    % DELAY-0
    SearchGene_0_C1 = strcmp(GeneSearchInC1C2,Delay0_Cell1_FINAL(:,2));
    LocOfGene_0_C1 = find(SearchGene_0_C1 == 1);
    SearchGene_0_C2 = strcmp(GeneSearchInC1C2,Delay0_Cell2_FINAL(:,2));
    LocOfGene_0_C2 = find(SearchGene_0_C2 == 1);
    FitScore_0_C1 = Delay0_Cell1_FINAL(LocOfGene_0_C1,1);
    FitScore_0_C2 = Delay0_Cell2_FINAL(LocOfGene_0_C2,1);

    % DELAY-20
    SearchGene_20_C1 = strcmp(GeneSearchInC1C2,Delay20_Cell1_FINAL(:,2));
    LocOfGene_20_C1 = find(SearchGene_20_C1 == 1);
    SearchGene_20_C2 = strcmp(GeneSearchInC1C2,Delay20_Cell2_FINAL(:,2));
    LocOfGene_20_C2 = find(SearchGene_20_C2 == 1);
    FitScore_20_C1 = Delay20_Cell1_FINAL(LocOfGene_20_C1,1);
    FitScore_20_C2 = Delay20_Cell2_FINAL(LocOfGene_20_C2,1);

    % DELAY-40
    SearchGene_40_C1 = strcmp(GeneSearchInC1C2,Delay40_Cell1_FINAL(:,2));
    LocOfGene_40_C1 = find(SearchGene_40_C1 == 1);
    SearchGene_40_C2 = strcmp(GeneSearchInC1C2,Delay40_Cell2_FINAL(:,2));
    LocOfGene_40_C2 = find(SearchGene_40_C2 == 1);
    FitScore_40_C1 = Delay40_Cell1_FINAL(LocOfGene_40_C1,1);
    FitScore_40_C2 = Delay40_Cell2_FINAL(LocOfGene_40_C2,1);

    % DELAY-60
    SearchGene_60_C1 = strcmp(GeneSearchInC1C2,Delay60_Cell1_FINAL(:,2));
    LocOfGene_60_C1 = find(SearchGene_60_C1 == 1);
    SearchGene_60_C2 = strcmp(GeneSearchInC1C2,Delay60_Cell2_FINAL(:,2));
    LocOfGene_60_C2 = find(SearchGene_60_C2 == 1);
    FitScore_60_C1 = Delay60_Cell1_FINAL(LocOfGene_60_C1,1);
    FitScore_60_C2 = Delay60_Cell2_FINAL(LocOfGene_60_C2,1);

    % DELAY-80
    SearchGene_80_C1 = strcmp(GeneSearchInC1C2,Delay80_Cell1_FINAL(:,2));
    LocOfGene_80_C1 = find(SearchGene_80_C1 == 1);
    SearchGene_80_C2 = strcmp(GeneSearchInC1C2,Delay80_Cell2_FINAL(:,2));
    LocOfGene_80_C2 = find(SearchGene_80_C2 == 1);
    FitScore_80_C1 = Delay80_Cell1_FINAL(LocOfGene_80_C1,1);
    FitScore_80_C2 = Delay80_Cell2_FINAL(LocOfGene_80_C2,1);

    % DELAY-100
    SearchGene_100_C1 = strcmp(GeneSearchInC1C2,Delay100_Cell1_FINAL(:,2));
    LocOfGene_100_C1 = find(SearchGene_100_C1 == 1);
    SearchGene_100_C2 = strcmp(GeneSearchInC1C2,Delay100_Cell2_FINAL(:,2));
    LocOfGene_100_C2 = find(SearchGene_100_C2 == 1);
    FitScore_100_C1 = Delay100_Cell1_FINAL(LocOfGene_100_C1,1);
    FitScore_100_C2 = Delay100_Cell2_FINAL(LocOfGene_100_C2,1);

    % COLLECT FITNESS SCORE
    Finess_C1 = [FitScore_0_C1 FitScore_20_C1 FitScore_40_C1...
                 FitScore_60_C1 FitScore_80_C1 FitScore_100_C1];

    Finess_C2 = [FitScore_0_C2 FitScore_20_C2 FitScore_40_C2...
                 FitScore_60_C2 FitScore_80_C2 FitScore_100_C2];
             
    clear Fitness_C1 Fitness_C2
    Fitness_C1 = [repmat('Fit=',size(InfoOnPlot_C1,1),1),...
                  char(num2str(round(cell2mat(Finess_C1(:,1)),1)))];   
    Fitness_C2 = [repmat('Fit=',size(InfoOnPlot_C2,1),1),...
                  char(num2str(round(cell2mat(Finess_C2(:,1)),1)))];
                          
    rf = 3;
    cf = 2;
else     
    Finess_C1 = Fitness_C1;
    Finess_C2 = Fitness_C2;
    rf = 2;
    cf = 2;
end     

%% ==========================         
%%  Plot figures and flip (inverse) if it is a 'REPRESSOR'
if FlipRepressor == 1
    ActRep_C1 = cellstr(ActRepInfo_C1);  
    ActRep_C2 = cellstr(ActRepInfo_C2);

    % To flip the repressor genes
    RepLocC1 = strcmp(', REP',ActRep_C1);   % Finding if the gene is repressor
    CheckIfAllRepC1 = all(RepLocC1 == 1);   % Checking if all probes are repressors in C1
    RepLocC2 = strcmp(', REP',ActRep_C2);
    CheckIfAllRepC2 = all(RepLocC2 == 1);

    if (CheckIfAllRepC1 == 1)  && (CheckIfAllRepC2 == 1)
        y1NORM = (y1NORM * (-1))+1; 
        y2NORM = (y2NORM * (-1))+1; 
        Flipped = ' {\color{red}FLIPPED}';
    else
        y1NORM = y1NORM;
        y2NORM = y2NORM;
        Flipped = '';
    end
else
    Flipped = '';
end



%% =============          PLOT FIGURES          =============
    % CELL-1
    if DisplayDelay == 1
        % To display max fit.score related delay information                 
        DelayInfo_C1 = [repmat(', ',size(InfoOnPlot_C1,1),1),...
                        repmat(char(strrep(MaxFitForDelay_C1,'-','')),size(InfoOnPlot_C1,1),1)];
        Legend_C1 = [Fitness_C1, FOXP3probe_C1, DelayInfo_C1, ProbeIDs_C1, PublicIDs_C1];
    else
        Legend_C1 = [Fitness_C1, FOXP3probe_C1, ProbeIDs_C1, PublicIDs_C1];%, CL_C1, ActRepInfo_C1];  
%     Legend_C1 = [Fitness_C1, FOXP3probe_C1, ProbeIDs_C1, CL_C1, ActRepInfo_C1];  
    end

    figure;
    subplot(rf,cf,1)
    % UNMORMALISED
    plot(time_points,GeneData_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
    hold on;
    plot(time_points,GeneData_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
    hold on;
    plot(time_points,DataToPlot_C1,'-*','Linewidth',1.8);
    maxInProbes_C1 = max([GeneData_C1(:);DataToPlot_C1(:)]);   
    xlim([0 360]); ylim([0 (maxInProbes_C1+maxInProbes_C1*(1/2))])
    xlabel('Time (min)'); ylabel('mRNA conc.'); 
    title(sprintf('%s', char(GeneSearchInC1C2), ' in Cell-1',... %num2str(cell2mat(InfoOnPlot_C1(1,8))),...
                  ' (',char(CodingNonCodFLAG),') ', Grade), 'FontSize',14);
    MyListForLegend_C1 = {'FOXP3 Probe-2', 'FOXP3 Probe-3'};
    for iii = 1:size(Legend_C1,1)
        MyListForLegend_C1{iii+2} = Legend_C1(iii,:);
    end
    legend(MyListForLegend_C1, 'FontWeight','bold', 'FontSize',11);
    
    % NORMALISED
    subplot(rf,cf,3)
    plot(time_points,GeneDataNORM_C1(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
    hold on;
    plot(time_points,GeneDataNORM_C1(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
    plot(time_points,y1NORM,'-*','Linewidth',1.8);
    xlim([0 360]); xlabel('Time (min)'); 
    ylabel('NORM. mRNA conc.');
    title(sprintf('%s', char(GeneSearchInC1C2), ' in Cell-1',... 
                  ' (NORMALISED)', Flipped),'FontSize',14); 
    
%%===================================================    
% CELL-2
    if DisplayDelay == 1
        % To display max fit.score related delay information                 
        DelayInfo_C2 = [repmat(', ',size(InfoOnPlot_C2,1),1),...
                        repmat(char(strrep(MaxFitForDelay_C2,'-','')),size(InfoOnPlot_C2,1),1)];
        Legend_C2 = [Fitness_C2, FOXP3probe_C2, DelayInfo_C2, ProbeIDs_C2, PublicIDs_C2];
    else
        Legend_C2 = [Fitness_C2, FOXP3probe_C2, ProbeIDs_C2, PublicIDs_C2];%, CL_C2, ActRepInfo_C2];  
%     Legend_C2 = [Fitness_C2, FOXP3probe_C2, ProbeIDs_C2, CL_C2, ActRepInfo_C2];  
    end   

    subplot(rf,cf,2)
    % UNMORMALISED
    plot(time_points,GeneData_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
    hold on;
    plot(time_points,GeneData_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
    hold on;
    plot(time_points,DataToPlot_C2,'-*','Linewidth',1.8);
    maxInProbes_C2 = max([GeneData_C2(:);DataToPlot_C2(:)]);   
    xlim([0 360]); ylim([0 (maxInProbes_C2+maxInProbes_C2*(1/2))])
    xlabel('Time (min)'); ylabel('mRNA conc.'); 
    title(sprintf('%s', char(GeneSearchInC1C2), ' in Cell-2',... 
                  ' (',char(CodingNonCodFLAG),') ', Grade), 'FontSize',14);
    MyListForLegend_C2 = {'FOXP3 Probe-2', 'FOXP3 Probe-3'};
    for jjj = 1:size(Legend_C2,1)
        MyListForLegend_C2{jjj+2} = Legend_C2(jjj,:);
    end
    legend(MyListForLegend_C2, 'FontWeight','bold', 'FontSize',11);

    % NORMALISED
    subplot(rf,cf,4)
    plot(time_points,GeneDataNORM_C2(1,:),'r--o','Color',[255 111 50]./255,'Linewidth',1.7);
    hold on;
    plot(time_points,GeneDataNORM_C2(2,:),'r--o','Color',[102 102 255]./255,'Linewidth',1.7);
    hold on;
    plot(time_points,y2NORM,'-*','Linewidth',1.8);
    xlim([0 360]); xlabel('Time (min)'); 
    ylabel('NORM. mRNA conc.');
    title(sprintf('%s', char(GeneSearchInC1C2), ' in Cell-2',... 
                  ' (NORMALISED)', Flipped),'FontSize',14);

              
%%===================================================    
% HISTOGRAM

if DelayHistogram == 1
    hT_C1 = [];
    hT_C2 = [];
    MyListForLegend_C1 = [];
    for jj = 1:size(Finess_C1,1)   % no. of occurences of each gene
        subplot(rf,cf,5)
        ProbeNum_C1 = char(InfoOnPlot_C1(jj,8));
        DispVal_C1 = [char(repmat(ProbeNum_C1(1:2),size(Finess_C1,2),1)) ...
                  char(repmat(' / ',size(Finess_C1,2),1))...
                  char(num2str(round(cell2mat(Finess_C1(jj,:)'),1)))];   
        b1 = bar(Delay,cell2mat(Finess_C1(jj,:)));
        b1.FaceAlpha = 0.4; % b1c = get(b1, 'Children');    
        hold on; grid on;
        xlabel('Delay (min)'); ylabel('Fitness Score')
        hT_C1=[hT_C1,text(b1(:).XData+b1(:).XOffset,b1(:).YData,DispVal_C1(1:6,:), ...
              'VerticalAlignment','bottom','horizontalalign','center','FontWeight','bold')];
    % %     labels = arrayfun(@(value) num2str(value,'%2.2f'),DispVal,'UniformOutput',false);
    % %     text(Delay,cell2mat(Finess_C1(jj,:)),labels,'HorizontalAlignment','center','VerticalAlignment','bottom')    
        title(['Cell-1, ' GeneSearchInC1C2{:}], 'FontWeight','bold', 'FontSize',13)
        %         legend([ab1{:}],Finess_C1(jj,:));

        subplot(rf,cf,6)
        ProbeNum_C2 = char(InfoOnPlot_C2(jj,8));
        DispVal_C2 = [char(repmat(ProbeNum_C2(1:2),size(Finess_C2,2),1)) ...
                  char(repmat(' / ',size(Finess_C2,2),1))...
                  char(num2str(round(cell2mat(Finess_C2(jj,:)'),1)))];   
        b2 = bar(Delay,cell2mat(Finess_C2(jj,:)),'r');
        b2.FaceAlpha = 0.4;
        hold on; grid on;
        xlabel('Delay (min)'); ylabel('Fitness Score','FontSize',12)
        hT_C2=[hT_C2,text(b2(:).XData+b2(:).XOffset,b2(:).YData,DispVal_C2(1:6,:), ...
              'VerticalAlignment','bottom','horizontalalign','center','FontWeight','bold')];
    % %     labels = arrayfun(@(value) num2str(value,'%2.2f'),cell2mat(Finess_C2(jj,:)),'UniformOutput',false);
    % %     text(Delay,cell2mat(Finess_C2(jj,:)),labels,'HorizontalAlignment','center','VerticalAlignment','bottom') 
        title(['Cell-2, ' GeneSearchInC1C2{:}], 'FontWeight','bold', 'FontSize',13)
    end
end


        % To Plot enlarged figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0.5 0.75 0.9]);





%************************       END       ******************************





  