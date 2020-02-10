
%%%% Script to filter out all time-series marked with "ABSENT" flag %%%%%

load('AFFIMETRIC_FLAGS_DataMatrixRawRecords.mat'); % probe IDs
load('AFFIMETRIC_FLAGS_DetectionCellArray.mat'); % Flags

ListOfAllFLAGS = [];
ListOfAllFLAGSSCORES = [];
for i = 1:numel(DataMatrixRawRecords.RowNames(:))
    %FlagLine = [];
        
    if  (i/1000-fix(i/1000)) == 0
        disp(i); 
    end
    
    Flag_Teff_Cell1 = 'NotABSENT';
    Flag_Teff_Cell2 = 'NotABSENT';
    Flag_Treg_Cell1 = 'NotABSENT';
    Flag_Treg_Cell2 = 'NotABSENT';
    
    %%%%% FLAG SCORE %%%%%
    FlagScore_Teff_Cell1 = 0;
    FlagScore_Teff_Cell2 = 0;
    FlagScore_Treg_Cell1 = 0;
    FlagScore_Treg_Cell2 = 0;
    %%%%%%%%%%%%%%%%%%%%%%
    
    NumberOfAbsentFLagsTeff1 = 0;
    NumberOfAbsentFLagsTeff2 = 0;
    NumberOfAbsentFLagsTreg1 = 0;
    NumberOfAbsentFLagsTreg2 = 0;
        
    for j = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19]
        
        MyCell = DetectionCellArray(i,j);
        if not(isempty(strfind(MyCell{1},'Absent'))) 
            NumberOfAbsentFLagsTeff1 = NumberOfAbsentFLagsTeff1 + 1;
        end
        
        MyCell = DetectionCellArray(i,19+j);
        if not(isempty(strfind(MyCell{1},'Absent')))
            NumberOfAbsentFLagsTeff2 = NumberOfAbsentFLagsTeff2 + 1;
        end
        
        MyCell = DetectionCellArray(i,2*19+18+j:2*19+18+j);
        if not(isempty(strfind(MyCell{1},'Absent')))
            NumberOfAbsentFLagsTreg2 = NumberOfAbsentFLagsTreg2 + 1;
        end
        
    end
    
    for k = 1:18
        MyCell = DetectionCellArray(i,2*19+k:2*19+k);
        if not(isempty(strfind(MyCell{1},'Absent')))
            NumberOfAbsentFLagsTreg1 = NumberOfAbsentFLagsTreg1 + 1;
        end
    end
    
    %%%%% FLAG SCORE %%%%%
    for j = 1:19
        
        MyCell = DetectionCellArray(i,j);
        if not(isempty(strfind(MyCell{1},'Absent'))) 
            FlagScore_Teff_Cell1 = FlagScore_Teff_Cell1 + 0;
        elseif not(isempty(strfind(MyCell{1},'Marginal'))) 
            FlagScore_Teff_Cell1 = FlagScore_Teff_Cell1 + 0.5;
        elseif not(isempty(strfind(MyCell{1},'Present'))) 
            FlagScore_Teff_Cell1 = FlagScore_Teff_Cell1 + 1.0;
        end
        
        MyCell = DetectionCellArray(i,19+j);
        if not(isempty(strfind(MyCell{1},'Absent'))) 
            FlagScore_Teff_Cell2 = FlagScore_Teff_Cell2 + 0;
        elseif not(isempty(strfind(MyCell{1},'Marginal'))) 
            FlagScore_Teff_Cell2 = FlagScore_Teff_Cell2 + 0.5;
        elseif not(isempty(strfind(MyCell{1},'Present'))) 
            FlagScore_Teff_Cell2 = FlagScore_Teff_Cell2 + 1.0;
        end
        
        MyCell = DetectionCellArray(i,2*19+18+j:2*19+18+j);
        if not(isempty(strfind(MyCell{1},'Absent'))) 
            FlagScore_Treg_Cell2 = FlagScore_Treg_Cell2 + 0;
        elseif not(isempty(strfind(MyCell{1},'Marginal'))) 
            FlagScore_Treg_Cell2 = FlagScore_Treg_Cell2 + 0.5;
        elseif not(isempty(strfind(MyCell{1},'Present'))) 
            FlagScore_Treg_Cell2 = FlagScore_Treg_Cell2 + 1.0;
        end
        
    end
    for k = 1:18
        MyCell = DetectionCellArray(i,2*19+k:2*19+k);
        if not(isempty(strfind(MyCell{1},'Absent'))) 
            FlagScore_Treg_Cell1 = FlagScore_Treg_Cell1 + 0;
        elseif not(isempty(strfind(MyCell{1},'Marginal'))) 
            FlagScore_Treg_Cell1 = FlagScore_Treg_Cell1 + 0.5;
        elseif not(isempty(strfind(MyCell{1},'Present'))) 
            FlagScore_Treg_Cell1 = FlagScore_Treg_Cell1 + 1.0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%
    
    if NumberOfAbsentFLagsTeff1 == 18 
        Flag_Teff_Cell1 = 'YesABSENT';
    end
    if NumberOfAbsentFLagsTeff2 == 19 
        Flag_Teff_Cell2 = 'YesABSENT';
    end
    if NumberOfAbsentFLagsTreg1 == 18 
        Flag_Treg_Cell1 = 'YesABSENT';
    end
    if NumberOfAbsentFLagsTreg2 == 19 
        Flag_Treg_Cell2 = 'YesABSENT';
    end
    
    FlagLine = {Flag_Teff_Cell1, Flag_Teff_Cell2, Flag_Treg_Cell1, Flag_Treg_Cell2};
    ListOfAllFLAGS = [ListOfAllFLAGS; FlagLine];
    
    %%%%% FLAG SCORE %%%%%
FlagSCORELine = [FlagScore_Teff_Cell1, FlagScore_Teff_Cell2, FlagScore_Treg_Cell1, FlagScore_Treg_Cell2];
    ListOfAllFLAGSSCORES = [ListOfAllFLAGSSCORES; FlagSCORELine];
    %%%%%%%%%%%%%%%%%%%%%%%
end

disp(ListOfAllFLAGSSCORES);
%disp(ListOfAllFLAGS);

load('GeneExpressionFengDataWithTimepointsNamesProbeIDs.mat');

CellOfGeneNamesProbIDsTimepointsRightOrder{2};
DataMatrixRawRecords.RowNames{63};
strcmp(CellOfGeneNamesProbIDsTimepointsRightOrder{2}, DataMatrixRawRecords.RowNames{63});

% parpool(4)
% %NewListOfAllFLAGSrightOrder = [];
% NewListOfAllFLAGSrightOrder = cell(numel(ListOfAllFLAGS(:,1)),1);
% parfor lineDATAFILE = 2:10%numel(CellOfGeneNamesProbIDsTimepointsRightOrder(:,1))
%     disp(lineDATAFILE)
%     for lineFLAGS = 1:numel(ListOfAllFLAGS(:,1))
%         if strcmp(CellOfGeneNamesProbIDsTimepointsRightOrder{lineDATAFILE}, DataMatrixRawRecords.RowNames{lineFLAGS})
%             disp(CellOfGeneNamesProbIDsTimepointsRightOrder(lineDATAFILE))
%             disp(DataMatrixRawRecords.RowNames(lineFLAGS))
%             disp(' ')
%             NewListOfAllFLAGSrightOrder{lineDATAFILE-1} = ListOfAllFLAGS(lineFLAGS,:);
%             %NewListOfAllFLAGSrightOrder = [NewListOfAllFLAGSrightOrder; ListOfAllFLAGS(lineFLAGS,:)];
%             break
%         end
%     end
%     %disp(lineDATAFILE);
%     %disp(NewListOfAllFLAGSrightOrder);
% end

% wO_firtst_line = CellOfGeneNamesProbIDsTimepointsRightOrder(:,2:54676)

count = 0;
NewListOfAllFLAGSrightOrder = cell(numel(ListOfAllFLAGSSCORES(:,1)),1);
for lineDATAFILE = 2:numel(CellOfGeneNamesProbIDsTimepointsRightOrder(:,1))
    if  (lineDATAFILE/1000-fix(lineDATAFILE/1000)) == 0
        disp(lineDATAFILE); 
    end
    count = count+1;
    %disp(CellOfGeneNamesProbIDsTimepointsRightOrder{count+1})
    %disp(DataMatrixRawRecords.RowNames(:))
    % 'scmp' search for each unique Probe ID of 'exprgcrma1' in 'ProbeGeneSym', 
    % and provides location of that matching ID in 'ProbeGeneSym'
    scmp(:,count) = strmatch(CellOfGeneNamesProbIDsTimepointsRightOrder{count+1}, DataMatrixRawRecords.RowNames(:));
    
end 
scmp = scmp';  
% 'id_gene' provides a cell array (54676 X 80) with information of Probe ID, gene symbol, gene name etc. along with data; 
% collected together and in the Probe ID sequence as that of 'exprgcrma1.mat'.
NewListOfAllFLAGSSCORErightOrder = ListOfAllFLAGSSCORES(scmp,:);

% disp(NewListOfAllFLAGSrightOrder);

save ListAFFIMETRICflags.mat 'NewListOfAllFLAGSrightOrder';

FileName = sprintf(['ListAFFIMETRICflags.mat']);
save(FileName, 'NewListOfAllFLAGSSCORErightOrder');

%disp(NewListOfAllFLAGSrightOrder)

