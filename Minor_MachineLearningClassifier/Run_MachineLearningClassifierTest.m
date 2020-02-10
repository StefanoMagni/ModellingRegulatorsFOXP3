%% CLASSIFIER

close all; clear; clc;

load('ResultsOfManualGeneSelection.mat')

% Converting struct to cell format
SELECTED_indices = struct2cell(selected_genes_indices);
MAYBE_indices = struct2cell(maybe_genes_indices);
REJECTED_indices = struct2cell(rejected_genes_indices);

% % Accessing each set of indices within each set
% SELECTED_set = SELECTED_indices(1);
% MAYBE_set = MAYBE_indices(1);
% REJECTED_set = REJECTED_indices(1);


% load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'Tcell1DataTable')
% data_set = Tcell1DataTable;
% load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell1')
% NamesTable = TableOfExtractedGeneNameAndIDs_cell1;
% load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'LjungBoxTestQandCL_Tcell1')
% LjungBoxValues = LjungBoxTestQandCL_Tcell1;
% % load('ExtractedCutData_Teff_MIGHT.mat', 'TableOfExtractedAFFIMETRICflagLines')



for ii = 4      % no of sets
    % jj = 1;    % no. of genes in a set

    if ii == 1
        load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'Tcell1DataTable')
        data_set = Tcell1DataTable;
        load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell1')
        NamesTable = TableOfExtractedGeneNameAndIDs_cell1;
        load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'LjungBoxTestQandCL_Tcell1')
        LjungBoxValues = LjungBoxTestQandCL_Tcell1;
        load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'FILTERED_AFFI_FLAGS_Tcell1')
        FilteredAFFIFlags = FILTERED_AFFI_FLAGS_Tcell1(:,1);
    elseif ii == 2
        load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'Tcell2DataTable')
        data_set = Tcell2DataTable;
        load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell2')
        NamesTable = TableOfExtractedGeneNameAndIDs_cell2;
        load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'LjungBoxTestQandCL_Tcell2')
        LjungBoxValues = LjungBoxTestQandCL_Tcell2;
        load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'FILTERED_AFFI_FLAGS_Tcell2')
        FilteredAFFIFlags = FILTERED_AFFI_FLAGS_Tcell2(:,2);
    elseif ii == 3
        load('FILTERED_ExtractedCutData_TregALLgen.mat', 'Tcell1DataTable')
        data_set = Tcell1DataTable;
        load('FILTERED_ExtractedCutData_TregALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell1')
        NamesTable = TableOfExtractedGeneNameAndIDs_cell1;
        load('FILTERED_ExtractedCutData_TregALLgen.mat', 'LjungBoxTestQandCL_Tcell1')
        LjungBoxValues = LjungBoxTestQandCL_Tcell1;
        load('FILTERED_ExtractedCutData_TregALLgen.mat', 'FILTERED_AFFI_FLAGS_Tcell2')
        FilteredAFFIFlags = FILTERED_AFFI_FLAGS_Tcell1(:,3);
    elseif ii == 4
        load('FILTERED_ExtractedCutData_TregALLgen.mat', 'Tcell2DataTable')
        data_set = Tcell2DataTable;
        load('FILTERED_ExtractedCutData_TregALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell2')
        NamesTable = TableOfExtractedGeneNameAndIDs_cell2;
        load('FILTERED_ExtractedCutData_TregALLgen.mat', 'LjungBoxTestQandCL_Tcell2')
        LjungBoxValues = LjungBoxTestQandCL_Tcell2;
        load('FILTERED_ExtractedCutData_TregALLgen.mat', 'FILTERED_AFFI_FLAGS_Tcell2')
        FilteredAFFIFlags = FILTERED_AFFI_FLAGS_Tcell2(:,4);
    end
    
    A = SELECTED_indices(ii);
    B = MAYBE_indices(ii);
    C = REJECTED_indices(ii);

    % Creating vectors of ones and zeros- of the same size as SELECTED,REJECTED genes
    SELECTED_ones = ones((size(A{1},1)),1);  % SELECTED_zeros = zeros((size(A{1},1)),1);     % SELECTED GENES
    % MAYBE_ones = ones((size(B{1},1)),1);  MAYBE_zeros = zeros((size(B{1},1)),1);    % MAYBE GENES
    REJECTED_zeros = zeros((size(C{1},1)),1); % REJECTED_ones = ones((size(C{1},1)),1);    % REJECTED GENES
    
    
    list_of_selected_genes_timeseries = data_set(SELECTED_indices{ii},:);
    list_of_selected_genes_Names = NamesTable(SELECTED_indices{ii},:);
    list_of_selected_genes_LjungBox = LjungBoxValues(SELECTED_indices{ii},:); %LjungBoxTestQandCL_Tcell1(SELECTED_indices{ii},:);
    list_of_selected_genes_AFFI_Flags = FilteredAFFIFlags(SELECTED_indices{ii},:);

    % list_of_maybe_genes_timeseries = data_set(MAYBE_indices{ii},:);
    % list_of_maybe_genes_Names = NamesTable(MAYBE_indices{ii},:);
    % list_of_maybe_genes_LjungBox = LjungBoxTestQandCL_Tcell1(MAYBE_indices{ii},:);
    % %list_of_maybe_genes_Names = NamesTable(MAYBE_indices{ii},:);

    list_of_rejected_genes_timeseries = data_set(REJECTED_indices{ii},:);
    list_of_rejected_genes_Names = NamesTable(REJECTED_indices{ii},:);
    list_of_rejected_genes_LjungBox = LjungBoxValues(REJECTED_indices{ii},:); %LjungBoxTestQandCL_Tcell1(REJECTED_indices{ii},:);
    list_of_rejected_genes_AFFI_Flags = FilteredAFFIFlags(REJECTED_indices{ii},:);
    
    combined_data_table = [list_of_selected_genes_timeseries; list_of_rejected_genes_timeseries];
    combined_response = [SELECTED_ones; REJECTED_zeros];
    combined_LjungBox = [list_of_selected_genes_LjungBox; list_of_rejected_genes_LjungBox];
    combined_AFFI_Flags = [list_of_selected_genes_AFFI_Flags; list_of_rejected_genes_AFFI_Flags];
    
end

% load('ExtractedCutData_Teff_MIGHT.mat', 'Tcell1DataTable')
% load('ExtractedCutData_Teff_MIGHT.mat', 'LjungBoxTestQandCL_Tcell1')
% load('ExtractedCutData_Teff_MIGHT.mat', 'TableOfExtractedAFFIMETRICflagLines')

timeseries_coeffofvariation = var(combined_data_table(:,:)')'./mean(combined_data_table(:,:)')'; % log(var(combined_data_table(:,:)')'); % 
timeseries_LjungBoxTestStatCL = combined_LjungBox(:,2);
timeseries_LjungBoxTestStatQ = combined_LjungBox(:,1);
% timeseries_AFFI_FLAGS = TableOfExtractedAFFIMETRICflagLines(:,1);

% count = 0;
% for i = 1:size(timeseries_coeffofvariation,1)
%     count = count +1;
%     timeseries_FoldIncrease6points(count,:) = max(Tcell1DataTable(count,1:6))/min(Tcell1DataTable(count,1:6));
% end
% figure;
% plot(1:19,Tcell1DataTable(1:189,:))

% combined_response = [1 1 1 1 1 1 1 0 0 1 0 1 1 1 0 0 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 1 1 1 0 0 0 0 0 1 1 1 0 1 1 1 0 1 1 1 1 1 1 1 1 0 1 1 0 0 0 0 1 1 1 0 0 0 0 0 1 1 0 0 1 0 0 1 0 1 0 1 1 0 1 0 1 1 1 0 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 0 1 0 0 0 1 1 1 1 0 0 0 1 0 0 0 0 0 1 1 0 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 0 0 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1]';

% count = 0;
% for i = 1:189
%     count = count + 1;
%     
%     % Compute fold increase in first 6 points (till t = 120min)
%     timeseries_FoldIncrease6points(count,:) = max(Tcell1DataTable(count,1:6))/min(Tcell1DataTable(count,1:6));
%     % Compute maximum value of and element in a column
%     [maxval(i), maxindx(i)] = max(Tcell1DataTable(i,:));   
%     % Compute average value of the time series response 
%     avg(i) = mean(Tcell1DataTable(i,1:19));
%     
% end
% maxval = maxval';
% avg = avg';


curvatureList = [];
for k = 1:size(combined_data_table,1)
    avg_1 = mean(combined_data_table(k,1:6))/mean(combined_data_table(k,1:19));
    avg_2 = mean(combined_data_table(k,7:14))/mean(combined_data_table(k,1:19));
    avg_3 = mean(combined_data_table(k,15:19))/mean(combined_data_table(k,1:19));

    x = [50 180 310];
    y = [avg_1 avg_2 avg_3];

    format long
    p = polyfit(x,y,2);
    curvature = abs(p(1)*1000000);
    
    curvatureList = [curvatureList, curvature];

%     figure;
%     plot(x,y,'*');
%     hold on;
%     t = [0:20:360];
%     plot(t,Tcell1DataTable(k,:)/mean(Tcell1DataTable(k,1:19)));
% 
%     xFit = 0:360;
%     yFit = polyval(p,xFit);
%     hold on;
%     plot(xFit,yFit);
%     hold off;

end

%     fold_change_avg = max(yFit)/min(yFit)

Combination = [timeseries_LjungBoxTestStatQ timeseries_LjungBoxTestStatCL ... 
    timeseries_coeffofvariation curvatureList' combined_AFFI_Flags combined_response]; % timeseries_AFFI_FLAGS

% USE SIMPLE TREE WITH Q AND log(VARIANCE) (79.4%)
% USE QUADRATIC SVM WITH CL, log(VARIANCE), AFFI FLAG SCORE (82%)

disp('DONE')
