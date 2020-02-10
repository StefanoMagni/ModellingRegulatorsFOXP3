

function PlotDynamicsBiographs(SelectedGeneData,Filtered_CutIntersection_ProbeN,...
         Tcell1DataTable,Tcell2DataTable,GeneNamesCell12, ProbeIDsCell1,...
         PublicIDsCell1,ProbeIDsCell2, PublicIDsCell2,GeneOfInterest,ProbeNumber,PlotBioGraph)
    Cl1_list = [];
    Cl2_list = [];
    load('CodingNonCodingList.mat');
    CodingNonCodSHORTlist = [];
    time_points = 0:20:360;

    for line = 1:size(Filtered_CutIntersection_ProbeN,1)

        % Extracting data of filtered genes from the main file to plot dynamic response
        y1 = Tcell1DataTable(cell2mat(Filtered_CutIntersection_ProbeN(line,1)),:);
        y2 = Tcell2DataTable(cell2mat(Filtered_CutIntersection_ProbeN(line,1+6)),:);

        % Ljung Box test (Compute Q and CL)
        NLags = 15; 
        Q1 = ComputeLjungBoxTestStatistics(y1, NLags);
        CL1 = chi2cdf(Q1,NLags);
        Cl1_list = [Cl1_list, CL1];
        Q2 = ComputeLjungBoxTestStatistics(y2, NLags);
        CL2 = chi2cdf(Q2,NLags);
        Cl2_list = [Cl2_list, CL2];

        % Find info on 'coding/non-coding/other'
        FindCodingNonCoding = strcmp(Filtered_CutIntersection_ProbeN(line,3),...
            CodingNonCodingList(:,2));
        locationOfFindResults = find(FindCodingNonCoding);
        CodingNonCodFLAGtemp = CodingNonCodingList(locationOfFindResults,:);
        % If no info is found on coding/non-coding then just put '-' instead
        if isempty(locationOfFindResults)
            CodingNonCodFLAGtemp = '-';
        end
        CodingNonCodFLAG = unique(CodingNonCodFLAGtemp(:,1));
        CodingNonCodSHORTlist = [CodingNonCodSHORTlist; CodingNonCodFLAG];

        hFig=figure;
% % %         WindowBegins = 100*ones(round(max([y1(:);y2(:)])));
% % %         WindowEnds = 240*ones(round(max([y1(:);y2(:)])));      
% % %         WindowLinesHeight = 1:round(max([y1(:);y2(:)]));
        
        str = sprintf([char(GeneNamesCell12(line,:)),'(',char(CodingNonCodFLAG),')']);

%         hh = suptitle([str ' - for ' num2str(GeneOfInterest)...
%          ' ProbeSet ID ' num2str(ProbeNumber) ' in Tregs']);
%          set(hh,'FontSize',20,'FontWeight','normal');

        subplot(2,1,1)%(3,1, 1 : 2)
        maxInCellsAllGenes = max([y1(:);y2(:)]); 
        y1NORM = (y1 - min(y1)) / (max(y1) - min(y1));
        y2NORM = (y2 - min(y2)) / (max(y2) - min(y2));
        SelectedGeneDataNORM1 = (SelectedGeneData(1,:) - min(SelectedGeneData(1,:))) / (max(SelectedGeneData(1,:)) - min(SelectedGeneData(1,:)));
        SelectedGeneDataNORM2 = (SelectedGeneData(2,:) - min(SelectedGeneData(2,:))) / (max(SelectedGeneData(2,:)) - min(SelectedGeneData(2,:)));
        plot(time_points,y1NORM,'b-o','Linewidth',2.5,'MarkerSize',10); 
        % plot(time_points,y1,'Linewidth',1.5); 
        xlim([0 360]); ylim([0 1]);%ylim([0 (maxInCellsAllGenes+(maxInCellsAllGenes*(1/4)))])
        hold on;
        %plot(time_points,y2,'Linewidth',1.5); xlim([0 360])
        plot(time_points,SelectedGeneDataNORM1,'b--*','Linewidth',2.5,'MarkerSize',10); xlim([0 360])
        hold on;
        ax = gca;
        ax.FontSize = 20;
% % %         plot(WindowBegins,WindowLinesHeight,'r--','Linewidth',0.9);
% % %         plot(WindowEnds,WindowLinesHeight,'r--','Linewidth',0.9);


        
        % Replacing '_' with '-' in Probe IDs string
        ProbeIDsCell1CurrentLine = strrep(ProbeIDsCell1(line,:),'_','-');
        ProbeIDsCell2CurrentLine = strrep(ProbeIDsCell2(line,:),'_','-');
        tempA = num2str(cell2mat(Filtered_CutIntersection_ProbeN(line,2)));
        tempB = num2str(cell2mat(Filtered_CutIntersection_ProbeN(line,2+6)));
        Fitness1 = ['\bf{Fitness=',tempA,'}'];
        Fitness2 = ['\bf{Fitness=',tempB,'}'];
        h=legend([char(Filtered_CutIntersection_ProbeN(line,3)),', ', ProbeIDsCell1CurrentLine,', '...
         PublicIDsCell1(line,:), char(10), 'LB-CL=' num2str(CL1), ', ',Fitness1],...
         [num2str(GeneOfInterest) ' ProbeSet ID ' num2str(ProbeNumber)],...
         'Location','northwest'); %  ' - Treg 1'
%         h=legend([char(GeneNamesCell12(line,:)),', ', ProbeIDsCell1CurrentLine, char(10) , PublicIDsCell1(line,:),...
%                   ', ', 'LB-CL=' num2str(CL1), ', ',Fitness1],...
%                    [num2str(GeneOfInterest) ' ProbeSet ID ' num2str(ProbeNumber)],'Location','northwest');
%         h = legend({['blue' char(10) 'line'],'red line'});        

        Green = [144/256,239/256,144/256];
        Yellow = [248/256,222/256,137/256];
        %Red = [255/256,160/256,122/256];

        if strcmp(strtrim(ProbeIDsCell1CurrentLine),...
                strtrim(ProbeIDsCell2CurrentLine)) == 1
            RGBcolor = Green;
%         elseif strcmp(PublicIDsCell1(line,:),PublicIDsCell2(line,:)) == 1
%             RGBcolor = Yellow;
        else
            RGBcolor = Yellow;
        end
        set(h,'color',RGBcolor);
        title('Treg - Cell 1'); ylabel('NORM. mRNA concentration')

        % FOXP3 probe-1 in cell 1-2 (for small subplot)
%         load(FileNameImportantGene)   
        subplot(2,1,2)%(3,1, 3)
        plot(time_points,y2NORM,'r-o','Linewidth',2.5,'MarkerSize',10);
        hold on;
        plot(time_points,SelectedGeneDataNORM2,'r--*','Linewidth',2.5,'MarkerSize',10);
        maxInCells = max(SelectedGeneData(:,:)'); 
        maxInProbes = max(maxInCells);
        xlim([0 360]); ylim([0 1]); %(maxInProbes+maxInProbes*(1/4))])
        hold on;
% % %         plot(WindowBegins,WindowLinesHeight,'r--','Linewidth',0.9);
% % %         plot(WindowEnds,WindowLinesHeight,'r--','Linewidth',0.9);
        MyLegend = legend([char(Filtered_CutIntersection_ProbeN(line,3)),', ', ProbeIDsCell2CurrentLine, PublicIDsCell1(line,:),...
                  char(10), 'LB-CL=' num2str(CL2), ', ',Fitness2],...
                   [num2str(GeneOfInterest) ' ProbeSet ID ' num2str(ProbeNumber)]);
%         set(MyLegend,'color',Green);
        Green = [144/256,239/256,144/256];
        Yellow = [248/256,222/256,137/256];
        %Red = [255/256,160/256,122/256];

        if strcmp(strtrim(ProbeIDsCell1CurrentLine),...
                strtrim(ProbeIDsCell2CurrentLine)) == 1
            RGBcolor = Green;
%         elseif strcmp(PublicIDsCell1(line,:),PublicIDsCell2(line,:)) == 1
%             RGBcolor = Yellow;
        else
            RGBcolor = Yellow;
        end
        set(MyLegend,'color',RGBcolor);
        title('Treg - Cell 2');
        xlabel('time (min)'); ylabel('NORM. mRNA concentration')
        ax = gca;
        ax.FontSize = 20;
        % save figures
        temp=[ num2str(GeneOfInterest) '_Probe' num2str(ProbeNumber) '_',...
               num2str(line),'.eps']; 
        % save figure with colored background (solved)
        hFig.InvertHardcopy='off';   
        hFig.Color='white'; 
        %saveas(gcf,temp,'epsc');
        hh = suptitle([str ' - for ' num2str(GeneOfInterest)...
                  ' ProbeSet ID ' num2str(ProbeNumber) ' in Tregs']);
        set(hh,'FontSize',20,'FontWeight','normal');
        % Enlarge figure to full screen.
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        print(temp,'-depsc');

    end

% % %     %% BIOGRAPH
% % %     for Cell = 1:2
% % %         if PlotBioGraph == 1
% % %             threshold = 0.0;
% % %             fitness = cell2mat(Filtered_CutIntersection_ProbeN(:,2+(Cell-1)*6));
% % %             mRNAs_names = Filtered_CutIntersection_ProbeN(:,3+(Cell-1)*6);
% % %             FigName = [ num2str(GeneOfInterest) '_Probe' num2str(ProbeNumber) '_Biograph_' num2str(Cell)];
% % %             dcgains = cell2mat(Filtered_CutIntersection_ProbeN(:,4+(Cell-1)*6));
% % %             mRNA_name_SingleElement = [num2str(GeneOfInterest) 'Probe' num2str(ProbeNumber)];
% % % 
% % %             plot_biographMODforVectors(fitness,threshold,mRNAs_names,FigName,...
% % %                                        dcgains,mRNA_name_SingleElement); 
% % %         end
% % %     end    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%% DEFINE FUNCTION TO COMPUTE LJUNG BOX TEST STATISTIC %%%%%%%%%%
    function [TestStatistics] = ComputeLjungBoxTestStatistics(DataPointsVector, NLags)
    Ndatapoins = numel(DataPointsVector);
    MyGeneAutocorrelation = autocorr(DataPointsVector); 
    TestStatistics = Ndatapoins * (Ndatapoins+2) * ...
        sum(MyGeneAutocorrelation(2:NLags+1).^2 ./ ...
        ((Ndatapoins-1):-1:(Ndatapoins-NLags)));
    end
end