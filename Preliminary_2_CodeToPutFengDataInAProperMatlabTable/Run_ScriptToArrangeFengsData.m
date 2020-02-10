% Author: Stefano Magni, created on 17/12/2016, contact: stefano.magni@uni.lu
% Script to arrange Feng's data in a matlab cell with proper alignement
% between names (from HG-U133_Plus_2), probe ids and numerical values (from
% en.mat), using expr_gcrma to get the correct order of lines
clear; close all; clc;
tic; % start measuring running time of the script

% Extract all the content of file HG-U133_Plus_2-annot-revised.csv 
% as text and put in a matlab cell variable  
run('AutomaticallyGeneratedScriptToconvertCSVtoMat.m');
disp('Here 1')

% Extract all the content of file expr_gcrma.csv 
% as text and put in a matlab cell variable
run('AutomaticallyGeneratedScriptToconvertCSVtoMat2.m');

% Save the content of the above created variables in a file for further use
TempFileName = sprintf(['TemporaryVariablesStorage' '.mat']);
save(TempFileName, 'CellOfGeneNamesAndProbeIDs', 'CellWithProbeIDsInRightOrder');

Temp = size(CellOfGeneNamesAndProbeIDs);
NumberOfLinesA = Temp(1);

Temp = size(CellWithProbeIDsInRightOrder);
NumberOfLinesB = Temp(1);

CellOfGeneNamesProbIDsTimepointsRightOrder = CellOfGeneNamesAndProbeIDs;
disp('Here 2')

% 1) for line in cell 1 and line number
UpperEnd = NumberOfLinesA;
for i = 2:UpperEnd
    % Get Probe Id
    CurrentProbID = CellOfGeneNamesAndProbeIDs(i,1);
    
    % 2) look for that probe number in cell 2 and record the line number
    for j = 1:NumberOfLinesB
        ProbIDinCellRightOrder = CellWithProbeIDsInRightOrder(j,1);
        if strcmp(CurrentProbID{1},ProbIDinCellRightOrder{1})
            % 3) append the line number at the end of the corresp line in cell 1
            CellOfGeneNamesProbIDsTimepointsRightOrder{i,6} = num2str(j);
        end
    end

    if rem(i,100) == 0
        disp(['PHASE I: ', num2str(i), ' lines scanned, ', num2str(UpperEnd-i), ' remaining.' ]);
    end
end
disp('Here 3')

% Save the content of the above created variable in a file for further use
TempFileName = sprintf(['TemporaryVariablesStorage2' '.mat']);
save(TempFileName, 'CellOfGeneNamesProbIDsTimepointsRightOrder');

disp([' ']);
load('en.mat')
EmptyLinesFound = 0;
for i = 2:UpperEnd
    % 4) copy the line corresp to that line numebr from en.mat
    if not(isempty(CellOfGeneNamesProbIDsTimepointsRightOrder{i,6}))
        CorrectLineNumber = CellOfGeneNamesProbIDsTimepointsRightOrder{i,6};
        CurrentTimepointsLine = en(str2double(CorrectLineNumber),:); 
        % 5) substitute it at the place of the line number in cell 1
        for k = 1:76
            CellOfGeneNamesProbIDsTimepointsRightOrder{i,5+k} = CurrentTimepointsLine(k);
        end
    else
        EmptyLinesFound = EmptyLinesFound + 1;
    end
    if rem(i,100) == 0
        disp(['PHASE II: ', num2str(i), ' lines scanned, ', num2str(UpperEnd-i), ' remaining.' ]);
    end
end
disp('Here 4')

disp([' ']);
disp(['The number of lines for which a matching has not been found is: ', num2str(EmptyLinesFound)]);

NInfoCols = 5;
Ntimepoints = 19;
for i = NInfoCols+1 : NInfoCols+Ntimepoints
   CellOfGeneNamesProbIDsTimepointsRightOrder{1,i} = sprintf([ ...
       'GeneExpression(t=', num2str(20*(i-(NInfoCols+1))),'mins) - Teff cell 1']);
end
for i = NInfoCols+1+Ntimepoints : NInfoCols+2*Ntimepoints
   CellOfGeneNamesProbIDsTimepointsRightOrder{1,i} = sprintf([ ...
       'GeneExpression(t=', num2str(20*(i-(NInfoCols+1+Ntimepoints))),'mins) - Teff cell 2']);
end
disp('Here 5')

for i = NInfoCols+1+2*Ntimepoints : NInfoCols+3*Ntimepoints
   CellOfGeneNamesProbIDsTimepointsRightOrder{1,i} = sprintf([ ...
       'GeneExpression(t=', num2str(20*(i-(NInfoCols+1+2*Ntimepoints))),'mins) - Treg cell 1']);
end
for i = NInfoCols+1+3*Ntimepoints : NInfoCols+4*Ntimepoints
   CellOfGeneNamesProbIDsTimepointsRightOrder{1,i} = sprintf([ ...
       'GeneExpression(t=', num2str(20*(i-(NInfoCols+1+3*Ntimepoints))),'mins) - Treg cell 2']);
end
disp('Here 6')

% Save the content of the above created variable in a file for further use
FileName = sprintf(['GeneExpressionFengDataWithTimepointsNamesProbeIDs' '.mat']);
save(FileName, 'CellOfGeneNamesProbIDsTimepointsRightOrder');

% Now cell 1 should look like 5 columns with GeneName, probe id etc, and
% then 76 columns of timepoints from en.mat
% Now the produced file is ready to be eaten with the next script.

toc; % display duration of time for which the script have been running

