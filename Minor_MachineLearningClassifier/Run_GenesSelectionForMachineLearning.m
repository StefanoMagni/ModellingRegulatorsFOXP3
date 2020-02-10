

%% TO TEST DIFFERENT FEATURES OF TCELL RESPONSE 


clc;
close all;

SetToBeUsed = {'Teff1', 'Teff2', 'Treg1', 'Treg2'};

for iii = 1:size(SetToBeUsed,2)

    % load data here
    if iii == 1
        %load('ExtractedCutData_Teff_KNOWN.mat', 'Tcell1DataTable')
        load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'Tcell1DataTable')
        data_set = Tcell1DataTable;
    elseif  iii == 2
        %load('ExtractedCutData_Teff_KNOWN.mat', 'Tcell2DataTable')
        load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'Tcell2DataTable')
        data_set = Tcell2DataTable;
    elseif  iii == 3
        %load('ExtractedCutData_Treg_KNOWN.mat', 'Tcell1DataTable')
        load('FILTERED_ExtractedCutData_TregALLgen.mat', 'Tcell1DataTable')
        data_set = Tcell1DataTable;
    elseif  iii == 4
        %load('ExtractedCutData_Treg_KNOWN.mat', 'Tcell2DataTable')
        load('FILTERED_ExtractedCutData_TregALLgen.mat', 'Tcell2DataTable')
        data_set = Tcell2DataTable;
    end
    
    t = 0:20:360;

    %% MOVING AVERAGE

    curvatureList = [];
    count = 0;
    for k = 1:size(data_set,1)
        count = count+1;

        avg_1(k) = mean(data_set(k,1:6))/mean(data_set(k,1:19));  % average of section 1
        avg_2(k) = mean(data_set(k,7:14))/mean(data_set(k,1:19)); % average of section 2
        avg_3(k) = mean(data_set(k,15:19))/mean(data_set(k,1:19));% average of section 3
        x = [50 180 310];
        y = [avg_1(k) avg_2(k) avg_3(k)];

        format long
        p = polyfit(x,y,2);           % fit polynomial of degree 2 
        curvature = abs(p(1)*1000000);
        curvatureList = [curvatureList, curvature];
    end


    %% PLOT AND SELECT GENES
    clear count x y 
    count= 0;
    for j = 1:60:size(data_set,1)
        count = count + 1;

        w = waitforbuttonpress;
        x = [50 180 310];
        y = [avg_1(j) avg_2(j) avg_3(j)];

        if w == 0
            plot(x,y*mean(data_set(j,1:19)),'*')                                 % plot floating average
            hold on;
            plot(t,data_set(j,:),'Linewidth',1.5)  % plot gene dynamics

            p = polyfit(x,y,2);                           % fit polynomial of degree 2 
            xFit = 1:360;
            yFit = polyval(p,xFit);
            hold on;
            plot(xFit,yFit*mean(data_set(j,1:19)),'Linewidth',1.5)               % plot polyfit curve

            fold_change_avg = max(yFit)/min(yFit);        % fold change average
            sprintf('Gene no. = %g \nFold change = %g',j,fold_change_avg) 

            a = input('Accept this gene (y/n/e)? ','s')     % user input to accept or reject gene
            close(figure(1));

            if strcmpi(a,'y')
               select_genes(count,:) = 1;   % selected genes
            elseif strcmpi(a,'e')
               select_genes(count,:) = 2;   % maybe
            else
               select_genes(count,:) = 0;   % rejected genes
            end
            decision_list = select_genes;

        else
            disp('End')
        end
    end

    SELECTED_GENES_index = find(select_genes(:,1) == 1);      % index of selected genes
    Number_of_genes_selected = size(SELECTED_GENES_index,1)   % no. of selected genes

    MAYBE_GENES_index = find(select_genes(:,1) == 2);         % index of 'maybe' genes
    Number_of_genes_maybe = size(MAYBE_GENES_index,1)      % no. of 'maybe' genes

    REJECTED_GENES_index = find((select_genes(:,1) == 0));    % index of rejected genes
    Number_of_genes_rejected = size(REJECTED_GENES_index,1)   % no. of rejected genes

    selected_gene_index = {'selected_Teff1', 'selected_Teff2', 'selected_Treg1', 'selected_Treg2'};
    maybe_gene_index = {'maybe_Teff1', 'maybe_Teff2', 'maybe_Treg1', 'maybe_Treg2'};
    rejected_gene_index = {'rejected_Teff1', 'rejected_Teff2', 'rejected_Treg1', 'rejected_Treg2'};
    
    selected_genes_indices.(selected_gene_index{iii}) = SELECTED_GENES_index;
    maybe_genes_indices.(maybe_gene_index{iii}) = MAYBE_GENES_index;
    rejected_genes_indices.(rejected_gene_index{iii}) = REJECTED_GENES_index;

end

save('ResultsOfManualGeneSelection.mat', ...
     'selected_genes_indices', 'maybe_genes_indices', 'rejected_genes_indices')




