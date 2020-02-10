% Author: Stefano Magni, stefano.magni@uni.lu, created on 01.03.2017

function GeneDataFiltering( NameOfFileWithNOTFilteredData, ApplyFengFilters, ApplyLjungBoxFilter, ...
    ApplyAFFIMETRICflagsFilter, LjungBoxTestCLTrashold, ThrasholdFengFILTERING, ...
    NfoldFILTERING, ThrasholdMeanFILTERING, MinOverallDifferenceFILTERING, ...
    FilterProbesWithMultipleGenes, DrawPlots, CellTypeAndGenesSubsetForPlot, ...
    ThrasholdFengFILTERINGleftTail, ApplyBonferroniCorrection, Silent, NLags)
%GENEDATAFILTERING Function to filter out the data on gene expression which
%do not fullfill criteria to differentiate them significantly from noise

tic; % start measuring running time of the script

% Load file containing NOT-filtered data
load(NameOfFileWithNOTFilteredData);

FILTERED_TableOfExtractedGeneNameAndIDs_cell1 = [];
FILTERED_TableOfExtractedGeneNameAndIDs_cell2 = [];

FILTERED_Tcell1DataTable = [];
FILTERED_Tcell2DataTable = [];

FILTERED_LjungBoxTestQandCL_Tcell1 = [];
FILTERED_LjungBoxTestQandCL_Tcell2 = [];

FILTERED_AFFI_FLAGS_Tcell1 = [];
FILTERED_AFFI_FLAGS_Tcell2 = [];

if TregORTeff == 'Teff'
    SHIFT = 0;
elseif TregORTeff == 'Treg'
    SHIFT = 2;
else
    disp('Error in Treg/Teff flag from filtering script');
end

% loop over cells
for CellIndex = 1:2
    kkkk=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if CellIndex == 1
        data = Tcell1DataTable;
        LjungBoxTestQandCL = LjungBoxTestQandCL_Tcell1;
    elseif CellIndex == 2
        data = Tcell2DataTable;
        LjungBoxTestQandCL = LjungBoxTestQandCL_Tcell2;
    end
    % FileName = sprintf(['Results_', TcellType{1}, '_', KnownMight{1}, '_Cell' num2str(CellIndex) '.mat']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over genes (i.e. lines)
    j = 0;
    ListOfFILTERmeOUTTAGSbeforeLjungBoxTest = [];
    for line = 1:1:numel(data(:,1))
        
        FILTERmeOUT_TAG = 0;

        OneGeneExpressions = data(line,:);
        [valMAX idxMAX] = max(OneGeneExpressions);
        [valMIN idxMIN] = min(OneGeneExpressions);
        if not(Silent)
            disp(line);
        end
        % 1) VERIFY IF FILTERING CONDITIONS ARE SATISFIED
        Table4AFFIMETRICflags = TableOfExtractedAFFIMETRICflagLines(line,:); % AFFI FLAG
        AFFIMETRICflag = Table4AFFIMETRICflags(CellIndex + SHIFT); % AFFI FLAG
        if (not(isempty(strfind(TableOfExtractedGeneNameAndIDs{line,4},'/')))) && (FilterProbesWithMultipleGenes == 1)
             % Filter out probes corresponding to more then one gene
             FILTERmeOUT_TAG = 1;
             kkkk = kkkk + 1;
             if not(Silent)
                 disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' num2str(line) ' was filtered out due to FILTER ///!']);
             end
        elseif ApplyAFFIMETRICflagsFilter == 1 && AFFIMETRICflag == 0 % strcmp(AFFIMETRICflag{1},'YesABSENT') % AFFI FLAG
            % Filter out all time-series gene expression data for which the
            % corresponding affimetric flag was 'Absent' in all the 19
            % datapoints.
            FILTERmeOUT_TAG = 1;
            if not(Silent)
                disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' ...
                num2str(line) ' was filtered out due to AFFIMETRIC FLAG because was ABSENT']);
            end
%         elseif valMAX < ThrasholdFengFILTERING && ApplyFengFilters
%             % Filter out everything smaller then given thrashold
%             FILTERmeOUT_TAG = 1;
%             disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' num2str(line) ' was filtered out due to FENG FILTER 1!']);
%         elseif valMAX / valMIN < NfoldFILTERING && ApplyFengFilters
%             % Filter out everything which is not changing N-fold
%             FILTERmeOUT_TAG = 1;
%             disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' num2str(line) ' was filtered out due to FENG FILTER 2!']);
        elseif (mean(OneGeneExpressions)<ThrasholdMeanFILTERING) || (valMAX<MinOverallDifferenceFILTERING) && ApplyFengFilters %-valMIN)
             % Filter out everything which is not changing N-fold
             FILTERmeOUT_TAG = 1;
             if not(Silent)
                disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' num2str(line) ' was filtered out due to FENG FILTER!']);
             end
%                 %&& (1-LjungBoxTestCLTrashold)/2.0 < LjungBoxTestStatisticsCL
%             % Filter out everything with a Q of the Ljung Box test
%             % corresponding to a CL lower then desired CL (CONFIDENCE LEVEL)
%             FILTERmeOUT_TAG = 1;
%             disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' ...
%                 num2str(line) ' was filtered out due to LjungBoxTest because Q = ' ...
%                 num2str(LjungBoxTestStatisticsQ) ' and Cumulative of X2 = ' ...
%                 num2str(LjungBoxTestStatisticsCL) ' while C.L. = ' ...
%                 num2str(LjungBoxTestCLTrashold) ' (2 tails test)']) ;
%         elseif ApplyLjungBoxFilter == 1 && LjungBoxTestStatisticsCL <= (1-LjungBoxTestCLTrashold)/2.0 ...
%                 && (not(isempty(strfind(TableOfExtractedGeneNameAndIDs{line,4},'/')))) 
%             % Filter out /// in left tail
%             FILTERmeOUT_TAG = 1;
%             disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' ...
%                 num2str(line) 'was filtered out as in left tail and having ///']) ;
%         elseif ApplyLjungBoxFilter == 1 && LjungBoxTestStatisticsCL <= (1-LjungBoxTestCLTrashold)/2.0 ...
%                && valMAX < ThrasholdFengFILTERINGleftTail
%             % Filter out < 500 in left tail
%             FILTERmeOUT_TAG = 1;
%             disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' ...
%                 num2str(line) 'was filtered out as in left tail and max < thrashold']) ;
%         elseif ApplyLjungBoxFilter == 1 && LjungBoxTestStatisticsCL <= (1-LjungBoxTestCLTrashold)/2.0 ...
%                && (not(isempty(strfind(TableOfExtractedGeneNameAndIDs{line,4},'---'))))
%              % Filter out probes corresponding to more then one gene
%              FILTERmeOUT_TAG = 1;
%              disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' num2str(line) ' was filtered out as in left tail and ---!!!']);
        else
            j = j + 1;
            if not(Silent)
                disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' num2str(line) ' was NOT FILTERED OUT before Ljung Box Test']);
            end
        end
        % disp('Here I am');
%         if not(isempty(strfind(TableOfExtractedGeneNameAndIDs{line,4},'---')))
%             TableOfExtractedGeneNameAndIDs{line,4}=TableOfExtractedGeneNameAndIDs{line,2}
%             disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' ...
%             num2str(line) 'was found to have --- and thus column 2 is used as name.']) ;
%         end
    ListOfFILTERmeOUTTAGSbeforeLjungBoxTest = [ListOfFILTERmeOUTTAGSbeforeLjungBoxTest, FILTERmeOUT_TAG]; 
    end
    jjj=j;    
    
    i=0;
    for line = 1:1:numel(data(:,1))
        
        FILTERmeOUT_TAG = ListOfFILTERmeOUTTAGSbeforeLjungBoxTest(line);
        
        OneGeneExpressions = data(line,:);
        
        format long
        LjungBoxTestStatisticsQ = LjungBoxTestQandCL(line,1);
        NDegOfFreedom = NLags;
        LjungBoxTestStatisticsCL = chi2cdf(LjungBoxTestStatisticsQ,NDegOfFreedom);%LjungBoxTestQandCL(line,2);         
        if not(Silent)
            disp(line);
        end
        
        if ApplyBonferroniCorrection == 0
            NsamplesBonferroniCorr = 1.;
        elseif ApplyBonferroniCorrection == 1
            NsamplesBonferroniCorr = jjj-kkkk;
        end
        
        if FILTERmeOUT_TAG == 0
            if ApplyLjungBoxFilter == 1 && ...
                LjungBoxTestStatisticsCL < 1-(1-LjungBoxTestCLTrashold)/(2.0*NsamplesBonferroniCorr) ...
                    && (1-LjungBoxTestCLTrashold)/(2.0*NsamplesBonferroniCorr) < LjungBoxTestStatisticsCL
                % Filter out everything with a Q of the Ljung Box test
                % corresponding to a CL lower then desired CL (CONFIDENCE LEVEL)
                FILTERmeOUT_TAG = 1;
                if not(Silent)
                    disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' ...
                        num2str(line) ' was filtered out due to LjungBoxTest because Q = ' ...
                        num2str(LjungBoxTestStatisticsQ) ' and Cumulative of X2 = ' ...
                        num2str(LjungBoxTestStatisticsCL) ' while C.L. = ' ...
                        num2str(LjungBoxTestCLTrashold) ' (1 tail test)']) ;
                end
            else
                i = i + 1;
                if not(Silent)
                    disp(['Gene ' TableOfExtractedGeneNameAndIDs{line,4} ' on line ' num2str(line) ' was NOT FILTERED OUT.']);
                end
            end
        end
        format short
        
        % 2) CUT ONLY FILTERED DATA
        if not(FILTERmeOUT_TAG)
            NewLine = [];
            NewLineBis = [];
            NewLineTris = [];
            NewLinePoker = [];
            if CellIndex == 1
                for j = 1:1:numel(TableOfExtractedGeneNameAndIDs(line,:))
                    NewLine = [NewLine, TableOfExtractedGeneNameAndIDs(line,j)];
                end
                FILTERED_TableOfExtractedGeneNameAndIDs_cell1 = ...
                    [FILTERED_TableOfExtractedGeneNameAndIDs_cell1; NewLine];
                for k = 1:1:numel(Tcell1DataTable(line,:))
                     NewLineBis = [NewLineBis, Tcell1DataTable(line,k)];
                end
                FILTERED_Tcell1DataTable = [FILTERED_Tcell1DataTable; NewLineBis];
                for y = 1:1:numel(LjungBoxTestQandCL(line,:))
                     NewLineTris = [NewLineTris, LjungBoxTestQandCL_Tcell1(line,y)];
                end
                FILTERED_LjungBoxTestQandCL_Tcell1 = [FILTERED_LjungBoxTestQandCL_Tcell1; NewLineTris];
                for j = 1:1:numel(TableOfExtractedAFFIMETRICflagLines(line,:))
                    NewLinePoker = [NewLinePoker, TableOfExtractedAFFIMETRICflagLines(line,j)];
                end
                FILTERED_AFFI_FLAGS_Tcell1 = ...
                    [FILTERED_AFFI_FLAGS_Tcell1; NewLinePoker];
            elseif CellIndex == 2
                for j = 1:1:numel(TableOfExtractedGeneNameAndIDs(line,:))
                    NewLine = [NewLine, TableOfExtractedGeneNameAndIDs(line,j)];
                end
                FILTERED_TableOfExtractedGeneNameAndIDs_cell2 = ...
                    [FILTERED_TableOfExtractedGeneNameAndIDs_cell2; NewLine];
                for k = 1:1:numel(Tcell2DataTable(line,:))
                     NewLineBis = [NewLineBis, Tcell2DataTable(line,k)];
                end
                FILTERED_Tcell2DataTable = [FILTERED_Tcell2DataTable; NewLineBis];
                for y = 1:1:numel(LjungBoxTestQandCL(line,:))
                     NewLineTris = [NewLineTris, LjungBoxTestQandCL_Tcell2(line,y)];
                end
                FILTERED_LjungBoxTestQandCL_Tcell2 = [FILTERED_LjungBoxTestQandCL_Tcell2; NewLineTris];
                for j = 1:1:numel(TableOfExtractedAFFIMETRICflagLines(line,:))
                    NewLinePoker = [NewLinePoker, TableOfExtractedAFFIMETRICflagLines(line,j)];
                end
                FILTERED_AFFI_FLAGS_Tcell2 = ...
                    [FILTERED_AFFI_FLAGS_Tcell2; NewLinePoker];
            else 
                disp('Problem in Cell index!!!');
            end
        end     
    end
    disp(['Initially, there were ' num2str(numel(data(:,1))-kkkk) ' genes'])
    disp(['After Affimetric Flag and Feng filter, ' num2str(jjj-kkkk) ' probe set IDs were left.'])
    disp(['Bonferroni correction = ' num2str(NsamplesBonferroniCorr)])
    disp(['After applying also Ljung Box Test filter, ' num2str(i-kkkk) ' probe set IDs were left.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE DATAPOINTS FOR ALL GENES AS FUNC OF t

ColGeneName = 4;
ListOfCellTypesAndGenesSubsets = {'Teff_KNOWN'; 'Teff_MIGHT'; 'Treg_KNOWN'; 'Treg_MIGHT'};
load('ResponsiveGenesNames.mat');

% Loop over the list of cell types and genes subsets
for index = 1:numel(ListOfCellTypesAndGenesSubsets)
    CellTypeAndGenesSubset=char(ListOfCellTypesAndGenesSubsets(index,1));
    disp(CellTypeAndGenesSubset)

    % 3) CHARGE THE DESIRED GENES NAMES LIST
    if CellTypeAndGenesSubset=='Teff_KNOWN'
        TregORTeff = 'Teff'; 
        ListOfGeneNames = Teff_KNOWNresp_GenesNames;
    elseif CellTypeAndGenesSubset=='Teff_MIGHT'
        TregORTeff = 'Teff';  
        ListOfGeneNames = Teff_MIGHTresp_GenesNames;
    elseif CellTypeAndGenesSubset=='Treg_KNOWN'
        TregORTeff = 'Treg';  
        ListOfGeneNames = Treg_KNOWNresp_GenesNames;
    elseif CellTypeAndGenesSubset=='Treg_MIGHT'
        TregORTeff = 'Treg';
        ListOfGeneNames = Treg_MIGHTresp_GenesNames;
    elseif or(CellTypeAndGenesSubset=='TregCustom',CellTypeAndGenesSubset=='TeffCustom')
        TregORTeff = CustomChoiceOfTregOrTeff; 
        ListOfGeneNames = CustomListOfGeneNames; 
    else
        disp('Error in genesDataExtractor 2)!');    
    end
end

if DrawPlots == 1
    TimePointsList = [0.]; % (minutes)
    for i = 1:(size(FILTERED_Tcell1DataTable,2)-1)
        TimePointsList = [ TimePointsList, TimePointsList(i) + 20.];
    end

    figure;
    % Make subplot 1
    subplot(2,1,1);
    for i = 1:size(FILTERED_Tcell1DataTable,1)
        plot(TimePointsList,FILTERED_Tcell1DataTable(i,:),'-o'); hold on;
    end
    title('Cell 1'); xlabel('time (min)'); ylabel('mRNA concentration (a.u.?)');

    % Make legend 1
    if size(FILTERED_TableOfExtractedGeneNameAndIDs_cell1,1) < 80
        ListOfLegendNames = [];
        for j = 1:size(FILTERED_TableOfExtractedGeneNameAndIDs_cell1,1)
            ListOfLegendNames = [ListOfLegendNames, FILTERED_TableOfExtractedGeneNameAndIDs_cell1(j,ColGeneName)];
        end
        legend(ListOfLegendNames);
        lh=findall(gcf,'tag','legend');
        set(lh,'location','northeastoutside'); 
    end

    % Make subplot 2
    subplot(2,1,2)
    for i = 1:size(FILTERED_Tcell2DataTable,1)
        plot(TimePointsList,FILTERED_Tcell2DataTable(i,:),'-o'); hold on;
    end
    title('Cell 2'); xlabel('time (min)'); ylabel('mRNA concentration (a.u.?)');
    hold on;

    % Make legend 2
    if size(FILTERED_TableOfExtractedGeneNameAndIDs_cell2,1) < 80
        ListOfLegendNames = [];
        for j = 1:size(FILTERED_TableOfExtractedGeneNameAndIDs_cell2,1)
            ListOfLegendNames = [ListOfLegendNames, FILTERED_TableOfExtractedGeneNameAndIDs_cell2(j,ColGeneName)];
        end
        legend(ListOfLegendNames);
        lh=findall(gcf,'tag','legend');
        set(lh,'location','northeastoutside'); 
    end
    hold off;

    % Make general plot title
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Gene expression time-series for two ', TregORTeff, ' cells.'], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

    % Give a proper name to the figure and save it
    FigName = sprintf(['FILTERED_GenesExpression_' CellTypeAndGenesSubsetForPlot]);
    savefig(FigName);
elseif DrawPlots == 0
    disp('No Plots Drawn.')
else
    disp('Error in DrawPlot!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the filtered data in a new file
newNameTabCell1 = 'TableOfExtractedGeneNameAndIDs_cell1';
newNameTabCell2 = 'TableOfExtractedGeneNameAndIDs_cell2';
newNameCell1 = 'Tcell1DataTable';
newNameCell2 = 'Tcell2DataTable';
newNameLjungCell1 = 'LjungBoxTestQandCL_Tcell1';
newNameLjungCell2 = 'LjungBoxTestQandCL_Tcell2';

STab.(newNameTabCell1) = FILTERED_TableOfExtractedGeneNameAndIDs_cell1;
STab.(newNameTabCell2) = FILTERED_TableOfExtractedGeneNameAndIDs_cell2;

SCell1.(newNameCell1) = FILTERED_Tcell1DataTable;
SCell2.(newNameCell2) = FILTERED_Tcell2DataTable;

SLjungCell1.(newNameLjungCell1) = FILTERED_LjungBoxTestQandCL_Tcell1;
SLjungCell2.(newNameLjungCell2) = FILTERED_LjungBoxTestQandCL_Tcell2;

FileName = sprintf(['FILTERED_' NameOfFileWithNOTFilteredData]);

save(FileName, '-struct', 'STab');
save(FileName, '-struct', 'SCell1', '-append');
save(FileName, '-struct', 'SCell2', '-append');
save(FileName, '-struct', 'SLjungCell1', '-append');
save(FileName, '-struct', 'SLjungCell2', '-append');
save(FileName, 'ListOfGeneMatches', '-append');
save(FileName, 'TregORTeff', '-append');
save(FileName, 'INFOonTableOfExctractedGeneNameAndIDs', '-append');
save(FileName, 'INFOTcell1Table', '-append');
save(FileName, 'INFOTcell2Table', '-append');
save(FileName, 'FILTERED_AFFI_FLAGS_Tcell1', '-append');
save(FileName, 'FILTERED_AFFI_FLAGS_Tcell2', '-append');

%save(NameOfFileWithNOTFilteredData, 'TableOfExtractedGeneNameAndIDs', '-append');

toc;

end

