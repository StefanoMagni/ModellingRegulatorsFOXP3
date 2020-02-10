

function PlotDynamicsBiographsUNION(FileNameGeneOfInterestDynamics, Filtered_CutIntersection_ProbeN,...
         Tcell1DataTable,Tcell2DataTable,GeneNamesCell12, ProbeIDsCell1,...
         PublicIDsCell1,ProbeIDsCell2, PublicIDsCell2,GeneOfInterest,ProbeNumberCutOrdered,PlotBioGraph,CellType)
    Cl1_list = [];
    Cl2_list = [];
    load('CodingNonCodingList.mat');
    CodingNonCodSHORTlist = [];
    time_points = 0:20:360;

    for line = 1:size(Filtered_CutIntersection_ProbeN,1)

        ProbeNumberGeneOfInterestCell1 = ProbeNumberCutOrdered(line,1);
        ProbeNumberGeneOfInterestCell2 = ProbeNumberCutOrdered(line,2);
        
        if ProbeNumberGeneOfInterestCell1 == ProbeNumberGeneOfInterestCell2
            disp(['Skip line' num2str(line) ' since Prob-Cell1 = Probe-Cell2']);
        else
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
            subplot(3,1, 1 : 2)
            maxInCellsAllGenes = max([y1(:);y2(:)]); 
            plot(time_points,y1,'Linewidth',1.5); 
            xlim([0 360]); ylim([0 (maxInCellsAllGenes+(maxInCellsAllGenes*(1/4)))])
            hold on;
            plot(time_points,y2,'Linewidth',1.5); xlim([0 360])
            hold on;

            % Replacing '_' with '-' in Probe IDs string
            ProbeIDsCell1CurrentLine = strrep(ProbeIDsCell1(line,:),'_','-');
            ProbeIDsCell2CurrentLine = strrep(ProbeIDsCell2(line,:),'_','-');        
            tempA = num2str(cell2mat(Filtered_CutIntersection_ProbeN(line,2)));
            tempB = num2str(cell2mat(Filtered_CutIntersection_ProbeN(line,2+6)));
            Fitness1 = ['\bf{Fitness=',tempA,'}'];
            Fitness2 = ['\bf{Fitness=',tempB,'}'];
            h=legend(['Cell1, ', ProbeIDsCell1CurrentLine, ', ', PublicIDsCell1(line,:),...
                    ', LB-CL=' num2str(CL1), ', ',Fitness1],...
                   ['Cell2, ', ProbeIDsCell2CurrentLine, ', ', PublicIDsCell2(line,:),...
                    ', LB-CL=' num2str(CL2), ', ',Fitness2]);

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
            str = sprintf([char(GeneNamesCell12(line,:)),'(',char(CodingNonCodFLAG),')']);
            title(str); ylabel('mRNA concentration')

                %SelectedGeneName = GeneOfInterestProbesInCell12{IndexGeneName};
            VariableNameToLoadCell1 = [GeneOfInterest '_Probe' ...
                num2str(ProbeNumberGeneOfInterestCell1) '_' CellType '12']; 
            SelectedGeneDataCell1temp = cell2mat(struct2cell(...
                load(FileNameGeneOfInterestDynamics, VariableNameToLoadCell1)));
            SelectedGeneDataCell1 = SelectedGeneDataCell1temp(1,:);
            
            VariableNameToLoadCell2 = [GeneOfInterest '_Probe' ...
                num2str(ProbeNumberGeneOfInterestCell2) '_' CellType '12']; 
            SelectedGeneDataCell2temp = cell2mat(struct2cell(...
                load(FileNameGeneOfInterestDynamics, VariableNameToLoadCell2)));
            SelectedGeneDataCell2 = SelectedGeneDataCell2temp(2,:);

            % FOXP3 probe-1 in cell 1-2 (for small subplot)
            subplot(3,1, 3)
            plot(time_points,SelectedGeneDataCell1,'Linewidth',0.7); hold on;
            plot(time_points,SelectedGeneDataCell2,'Linewidth',0.7); hold off;
            maxInCell1 = max(SelectedGeneDataCell1(:,:)'); 
            maxInCell2 = max(SelectedGeneDataCell2(:,:)'); 
            maxInProbes = max(maxInCell1,maxInCell2);
            xlim([0 360]); ylim([0 (maxInProbes+maxInProbes*(1/4))])
            hold on;
            MyLegend = legend(['Cell-1, Probe ' num2str(ProbeNumberGeneOfInterestCell1)],...
                ['Cell-2, Probe ' num2str(ProbeNumberGeneOfInterestCell2)],...
                'Location','northwest');
            set(MyLegend,'color',Yellow);
            title([num2str(GeneOfInterest) ' in ' CellType 's']);
            xlabel('time (min)'); ylabel('mRNA conc.')
            % save figures
            temp=[ num2str(GeneOfInterest) '_ProbesUNION_' num2str(line) '.eps'];
            hFig.InvertHardcopy='off';   
            hFig.Color='white'; 
            print(temp,'-depsc');
            % saveas(gcf,temp,'epsc');
        end
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