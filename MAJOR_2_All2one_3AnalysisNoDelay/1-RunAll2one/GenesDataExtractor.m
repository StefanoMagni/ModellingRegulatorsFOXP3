% Author: Stefano Magni, created on 16/12/2016, contact: stefano.magni@uni.lu
% Script to extract desired portions of datafiles (gene expression data)
% according to provided lists of gene names, and to data structure.


function [] = GenesDataExtractor(TakeGeneNamesList, ListOfGeneNames, TregTeff, DrawPlots)
% JUST TYPE THE FOLLOWING LINE IN THE COMMAND WINDOW TO EXTRACT THE 4 LISTS OF FENG'S DATA:
% GenesDataExtractor(0, {}, '', 1)

%tic; % start measuring running time of the script

load('ListAFFIMETRICflags.mat', 'NewListOfAllFLAGSSCORErightOrder')

%%%%%%%% ONLY NEXT LINES WERE TO BE MODIFIED TO CONTROL THE SCRIPT, BUT NOW IT IS A FUNCTION %%%%%%%%
GeneNamesFromCustomList=TakeGeneNamesList; % Choose:
% =2 to use ALL genes!!!
% =1 to provide the list of desired gene names by hand, or
% =0 to get the 4 lists contained in the file ResponsiveGenesNames
% SPECIFY HERE THE CUSTOM LIST OF GENE NAMES TO BE EXCTRACTED AND THE DESIRED TYPE OF CELL
CustomListOfGeneNames = ListOfGeneNames;% {'GeneName1'; 'GeneName2'; 'GeneName3'; 'GeneName11'; 'GeneWeirdName'};
CustomChoiceOfTregOrTeff = TregTeff; % 'Teff'; % Choose among 'Treg' or 'Teff'
GenerateFakeGeneExpressionData = 0; % 0 to use Feng's real data, 1 to create fake data for test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET PARAMETERS DESCRIBING THE STUCTURE OF THE INPUT DATA FILE
NameOfInputDataFile = 'GeneExpressionFengDataWithTimepointsNamesProbeIDs.mat';
Ninfo = 5; % number of columns containing informations instead of timepoints, at the beginning of file
ColGname = 4; % Specify in which column the Gene Name is
NofTimepoints = 19; % number of time-points for each gene
TotNofCells = 4;
m = Ninfo + TotNofCells*NofTimepoints; % number of columns of the datatfile
n = 54675+1; % 100+1 number of lines of the datatfile
% ENSURE NUMBER OF LINES FOR DATA FILE IS 54675+1 WHEN USING REAL FENG's DATA!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) LOAD FILE CONTAINING THE DATA WE WANT TO EXTRACT
Data = cell(n,m);
if GenerateFakeGeneExpressionData == 0
    load(NameOfInputDataFile);
    Data = CellOfGeneNamesProbIDsTimepointsRightOrder; % Copy Real Data in a temp variable
%%%%%%%%%%%%%%%%%%%% OR CREATE FAKE DATA FOR TEST, AVOID WHEN USING REAL DATA %%%%%%%%%%%%%%
elseif GenerateFakeGeneExpressionData == 1
    for l=1:m                                  
        Data{1,l} = ['Description' num2str(l)];                    
    end  
    for i = 2:n                                        
        for k=1:Ninfo                            
            Data{i,k} = ['SomeInfo' num2str(i-1)];  
        end                                          
        Data{i,ColGname} = ['GeneName' num2str(i-1)];  
        for j = Ninfo+1:m                                  
            Data{i,j} = rand + i - 1;                    
        end                                          
    end     
    Data{14,ColGname} = ['GeneName3  //// GeneName7'];  
    Data{15,ColGname} = ['GeneName55  /// GeneName3'];  
    Data{16,ColGname} = ['GeneName99  /// GeneName35'];  
    Data{17,ColGname} = ['GeneName11'];  
    Data{18,ColGname} = ['GeneName11'];  
    Data{19,ColGname} = ['GeneName234'];  
    save('FakeDataFileForTest.mat','Data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    disp('Error with GenerateFakeGeneExpressionData');
end

if GeneNamesFromCustomList == 1
    ListOfCellTypesAndGenesSubsets = {[CustomChoiceOfTregOrTeff 'Custom']};
elseif GeneNamesFromCustomList == 0
    % Prepare list with 4 lists of names provided by Feng, to loop over
    ListOfCellTypesAndGenesSubsets = {'Teff_KNOWN'; 'Teff_MIGHT'; 'Treg_KNOWN'; 'Treg_MIGHT'};
    % 2) LOAD LISTS OF GENE NAMES WE WANT TO SELECT
    load('ResponsiveGenesNames.mat');
elseif GeneNamesFromCustomList == 2
    ListOfCellTypesAndGenesSubsets = {[CustomChoiceOfTregOrTeff 'ALLgen']};
else
    disp('PROBLEM at beginning of GenesDataExtractor!');
end

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
    elseif or(CellTypeAndGenesSubset=='TregALLgen',CellTypeAndGenesSubset=='TeffALLgen')
        TregORTeff = CustomChoiceOfTregOrTeff; 
        ListOfGeneNames = CellOfGeneNamesProbIDsTimepointsRightOrder(:,4);
    else
        disp('Error in genesDataExtractor 2)!');    
    end
    
    TableOfExtractedGeneNameAndIDs = [];
    TableOfExtractedDataLines = [];
    ListOfMatchingLines = [];
    
    TableOfExtractedAFFIMETRICflagLines = [];  % AFFI FLAG
    load('ListAFFIMETRICflags.mat');  % AFFI FLAG
    
    %%%%%%%%% NOW THIS WE SHOULD SKIP IF YOU WANT ALL GENES %%%%%%%%%%
    if GeneNamesFromCustomList ~= 2
        
        % 4) SEARCH MACTHING GENE NAMES INTO THE DATA FILE
        for j = 1:numel(ListOfGeneNames) % loop over all gene names
            ExtractedGeneNameAndID = [];
            ExtractedDataLine = []; 
            NumberOfLinesFoundToMatchGeneName = 0;
            for i = 2:n % loop over all lines of datafile BUT first (the one with info)
                DesiredGeneName = char(ListOfGeneNames(j)); 
                GeneNamesInCurrentFileLine = Data{i,ColGname};
                % search DesiredGeneName into GeneNamesInCurrentFileLine
                IndexOfFoundName = strfind(GeneNamesInCurrentFileLine,DesiredGeneName);
                if not(isempty(IndexOfFoundName)) % Do only if the desired gene name have been found
                    % Before proceeding further, verify if the found line really
                    % contains the gene name we searched for, or a longer one 
                    % This to avoid e.g. to mistakenly take GeneName34 when GeneName3 is wanted.
                    % Also to avoid e.g. to mistakenly take GeneName34 when Name34 is wanted.
                    % PART I: VERIFY AFTER GENE NAME
                    % A] count number of letters of current GeneName
                    LenghtGeneName = numel(DesiredGeneName);
                    % B] identify position at which the Gene Name has been found
                    % and advance untill end of desired Gene Name
                    FirstPositionAfterName = IndexOfFoundName + LenghtGeneName;
                    EXCTRACTLINE = 0;
                    if numel(GeneNamesInCurrentFileLine)==FirstPositionAfterName-1
                    % C] if this is last character, KEEP LINE
                        EXCTRACTLINE = 1;
                    else
                        if numel(FirstPositionAfterName)>1
                            EXCTRACTLINE = 1;
                        elseif GeneNamesInCurrentFileLine(FirstPositionAfterName)==' '
                        % E] else, if this is followed by a space, KEEP LINE
                            EXCTRACTLINE = 1;
                        else
                        % F] else, if this is followed by another character, DO NOT keep line.
                            EXCTRACTLINE = 0;
                        end
                    end
                    if EXCTRACTLINE==1
                        % PART II: VERIFY BEFORE GENE NAME
                        % B] identify position at which the Gene Name has been found
                        % and go back untill beginning of desired Gene Name
                        LastPositionBeforeName = IndexOfFoundName - 1;
                        if LastPositionBeforeName == 0
                        % C] if this is first character, KEEP LINE
                            EXCTRACTLINE = 1;
                        else
                            if numel(FirstPositionAfterName)>1
                            EXCTRACTLINE = 1;
                            elseif GeneNamesInCurrentFileLine(LastPositionBeforeName)==' '
                            % E] else, if this is preceeded by a space, KEEP LINE
                                EXCTRACTLINE = 1;
                            else
                            % F] else, if this is preceeded by another character, DO NOT keep line.
                                EXCTRACTLINE = 0;
                            end
                        end
                    end

                    if EXCTRACTLINE
                        % 5) WHEN CORRECTLY MATCHING GENE NAME, THEN EXTRACT DATA LINE 
                        ExtractedGeneNameAndID = {Data{i,1:Ninfo}};
                        for k = 1:(m-Ninfo)
                            ExtractedDataLine = [ExtractedDataLine, Data{i,k+Ninfo}];
                        end
                        NumberOfLinesFoundToMatchGeneName = NumberOfLinesFoundToMatchGeneName + 1;

                        TableOfExtractedGeneNameAndIDs = [TableOfExtractedGeneNameAndIDs; ExtractedGeneNameAndID];
                        TableOfExtractedDataLines = [TableOfExtractedDataLines; ExtractedDataLine];
                        TableOfExtractedAFFIMETRICflagLines = [TableOfExtractedAFFIMETRICflagLines; NewListOfAllFLAGSSCORErightOrder(i-1,:)]; % AFFI FLAG
                        ExtractedGeneNameAndID = [];
                        ExtractedDataLine = []; 
                    end
                end
            end
            ListOfMatchingLines = [ListOfMatchingLines, NumberOfLinesFoundToMatchGeneName]; 
        end

        disp(TableOfExtractedAFFIMETRICflagLines); % AFFI FLAG
        
        % 6) VERIFY GENES MATCHING BETWEEN LIST OF DESIRED NAMES AND DATA FILE
        for l = 1:numel(ListOfMatchingLines)
            if ListOfMatchingLines(l) ~= 1
                if ListOfMatchingLines(l) > 1
                    disp(sprintf(['\nWarning in matching: ', num2str(ListOfMatchingLines(l)), ...
                        ' lines matching with gene name ', strjoin(ListOfGeneNames(l)) ' have been found!\n']));
                elseif ListOfMatchingLines(l) == 0
                    disp(sprintf(['\nPROBLEM in matching: ', num2str(ListOfMatchingLines(l)), ...
                        ' lines matching with gene name ', strjoin(ListOfGeneNames(l)) ' have been found!\n']));
                end
            end
        end
        ListOfGeneMatches = ListOfGeneNames;
        for h = 1:numel(ListOfGeneMatches)
            ListOfGeneMatches{h,2} = ListOfMatchingLines(h);
        end
    elseif GeneNamesFromCustomList == 2
        disp('Start filling TableOfExtractedDataLines')
        ExtractedGeneNameAndID = [];
        ExtractedDataLine = []; 
        %TableOfExtractedGeneNameAndIDs = [TableOfExtractedGeneNameAndIDs; ExtractedGeneNameAndID];
        TableOfExtractedDataLines = [];
        tic
        for i = 2:1:n
            ExtractedGeneNameAndID = {Data{i,1:Ninfo}};
            for k = 1:(m-Ninfo)
                ExtractedDataLine = [ExtractedDataLine, Data{i,k+Ninfo}];
            end
            TableOfExtractedDataLines = [TableOfExtractedDataLines; ExtractedDataLine];
            TableOfExtractedGeneNameAndIDs = [TableOfExtractedGeneNameAndIDs; ExtractedGeneNameAndID];
            TableOfExtractedAFFIMETRICflagLines = [TableOfExtractedAFFIMETRICflagLines; NewListOfAllFLAGSSCORErightOrder(i-1,:)]; % AFFI FLAG
            ExtractedGeneNameAndID = [];
            ExtractedDataLine = []; 
            if mod(i,1000) == 0
                disp(i)
                toc
            end
        end
        ListOfGeneMatches = ['ALL genes included'];
        disp(TableOfExtractedAFFIMETRICflagLines); % AFFI FLAG
    end
    
    %%%%%%%%% HERE WERE YOU SHOULD START IF YOU WANT ALL GENES %%%%%%%%%%
    disp('7777')
    % 7) CUT FOR Teff1/Teff2/Treg1/Treg2
    tic
    if TregORTeff == 'Teff'
        Tcell1DataTable = [];
        Tcell2DataTable = [];
        for i = 1:size(TableOfExtractedDataLines,1)
            for j = (1) : (NofTimepoints)   
                Tcell1DataTable(i,j) =  TableOfExtractedDataLines(i,j);
            end    
            for j = (NofTimepoints+1) : (2*NofTimepoints)
                Tcell2DataTable(i,j-NofTimepoints) =  TableOfExtractedDataLines(i,j);
            end
            if mod(i,1000) == 0
                disp('Now splitting Teff tables, line: '); disp(i);
                toc
            end
        end
    elseif TregORTeff == 'Treg'
        Tcell1DataTable = [];
        Tcell2DataTable = [];
        for i = 1:size(TableOfExtractedDataLines,1)
            for j = (2*NofTimepoints+1) : (3*NofTimepoints)   
                Tcell1DataTable(i,j-2*NofTimepoints) =  TableOfExtractedDataLines(i,j);
            end    
            for j = (3*NofTimepoints+1) : (4*NofTimepoints)
                Tcell2DataTable(i,j-3*NofTimepoints) =  TableOfExtractedDataLines(i,j);
            end
            if mod(i,1000) == 0
                disp('Now spplitting Treg tables, line: '); disp(i);
            end
        end
    else
        disp('The List Of Genes Names has an anknown name and then we cannot cut to separe treg and teff.')
    end

    % 8) Exctract info about column content, to be saved separately
    INFOonTableOfExctractedGeneNameAndIDs = cell(1,Ninfo);
    for i = 1:Ninfo
        INFOonTableOfExctractedGeneNameAndIDs{1,i} = Data{1,i};
    end
    InfoTcell1Table = cell(1,NofTimepoints);
    InfoTcell2Table = cell(1,NofTimepoints);
    for k = 1:NofTimepoints
        if TregORTeff == 'Teff'
            INFOTcell1Table{1,k} =  Data{1,k+Ninfo};
            INFOTcell2Table{1,k} =  Data{1,k+Ninfo+NofTimepoints};
        elseif TregORTeff == 'Treg'
            INFOTcell1Table{1,k} =  Data{1,k+Ninfo+2*NofTimepoints};
            INFOTcell2Table{1,k} =  Data{1,k+Ninfo+3*NofTimepoints}; 
        end
    end
    
    % 9) SAVE RESULTS FOR USAGE WITH ANOTHER SCRIPT
    FileName = sprintf(['ExtractedCutData_' CellTypeAndGenesSubset '.mat']);
    save(FileName, 'TableOfExtractedGeneNameAndIDs', 'Tcell1DataTable', ...
        'Tcell2DataTable', 'ListOfGeneMatches', 'TregORTeff', ...
        'INFOonTableOfExctractedGeneNameAndIDs', 'INFOTcell1Table', 'INFOTcell2Table', 'TableOfExtractedAFFIMETRICflagLines');

    % 10) PLOT THE DATAPOINTS FOR ALL GENES AS FUNC OF t
    if DrawPlots == 1
        TimePointsList = [0.]; % (minutes)
        for i = 1:(size(Tcell1DataTable,2)-1)
            TimePointsList = [ TimePointsList, TimePointsList(i) + 20.];
        end

        figure;
        % Make subplot 1
        subplot(2,1,1);
        for i = 1:size(Tcell1DataTable,1)
            plot(TimePointsList,Tcell1DataTable(i,:),'-o'); hold on;
        end
        title('Cell 1'); xlabel('time (min)'); ylabel('mRNA concentration (a.u.?)');

        % Make subplot 2
        subplot(2,1,2)
        for i = 1:size(Tcell2DataTable,1)
            plot(TimePointsList,Tcell2DataTable(i,:),'-o'); hold on;
        end
        title('Cell 2'); xlabel('time (min)'); ylabel('mRNA concentration (a.u.?)');
        hold on;

        % Make legend
        if size(TableOfExtractedGeneNameAndIDs,1) < 80
            ListOfLegendNames = [];
            for j = 1:size(TableOfExtractedGeneNameAndIDs,1)
                ListOfLegendNames = [ListOfLegendNames, TableOfExtractedGeneNameAndIDs(j,ColGname)];
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
        FigName = sprintf(['GenesExpression_' CellTypeAndGenesSubset]);
        savefig(FigName);
    elseif DrawPlots == 0
        disp('No Plots Drawn.')
    else
        disp('Error in DrawPlot!')
    end
end

% 11) CLEAR AND DISPLAY EXECUTION TIME OF THE SCRIPT
clear
%toc; % display duration of time for which the script has been running
end
