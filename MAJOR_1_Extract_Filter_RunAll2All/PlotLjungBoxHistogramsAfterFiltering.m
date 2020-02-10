
%disp('HEY!!!!')
MakeHistogramsOfLjungTestStatistics_FILTERED = 1;
MakeHistogramsOfLjungTestStatisticsALL_FILTERED = 1;

if MakeHistogramsOfLjungTestStatistics_FILTERED == 1
    
    ListOfListsToBePlot = {};
    for TcellType = {'Teff_KNOWN', 'Teff_MIGHT', 'Treg_KNOWN', 'Treg_MIGHT'}

        FileNameToBeLoaded = ['FILTERED_ExtractedCutData_' , TcellType{1}, '.mat'];
        load(FileNameToBeLoaded)
        for TableToBeUsed = {LjungBoxTestQandCL_Tcell1, LjungBoxTestQandCL_Tcell2}
            %disp(TableToBeUsed{1});
            %disp(TableToBeUsed{1}(:,1));
            ListOfListsToBePlot{end+1} = transpose(TableToBeUsed{1}(:,1));
        end

    end

    % Make Histograms of the distribution of the Test Statistics for all
    % Treg/Teff x cell1/cell2 x KNOWN/MIGHT
    figure

    subplot(2,2,1)
    hist(ListOfListsToBePlot{3},20);
    h = findobj(gca,'Type','patch');
    hold on
    hist(ListOfListsToBePlot{1},20);
    % Color the 2 histograms differently
    set(findobj(gca,'Type','patch'),'FaceColor','r','EdgeColor','w');
    set(h,'FaceColor','b','EdgeColor','w');
    title('Teff - Cell 1'); xlabel('Test statistics Q (Ljung-Box test)'); ylabel('Number of genes per bin');
    % Make legend
    ListOfLegendNames = {'FILTERED MIGHT BE responsive genes'; 'KNOWN responsive genes'};
    legend(ListOfLegendNames);
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeast'); 
    hold off

    subplot(2,2,2)
    hist(ListOfListsToBePlot{7},20);
    h = findobj(gca,'Type','patch');
    hold on
    hist(ListOfListsToBePlot{5},20);
    % Color the 2 histograms differently
    set(findobj(gca,'Type','patch'),'FaceColor','r','EdgeColor','w');
    set(h,'FaceColor','b','EdgeColor','w');
    title('Treg - Cell 1'); xlabel('Test statistics Q (Ljung-Box test)'); ylabel('Number of genes per bin');
    % Make legend
    ListOfLegendNames = {'MIGHT BE responsive genes'; 'KNOWN responsive genes'};
    legend(ListOfLegendNames);
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeast'); 
    hold off

    subplot(2,2,3)
    hist(ListOfListsToBePlot{4},20);
    h = findobj(gca,'Type','patch');
    hold on
    hist(ListOfListsToBePlot{2},20);
    % Color the 2 histograms differently
    set(findobj(gca,'Type','patch'),'FaceColor','r','EdgeColor','w');
    set(h,'FaceColor','b','EdgeColor','w');
    title('Teff - Cell 2'); xlabel('Test statistics Q (Ljung-Box test)'); ylabel('Number of genes per bin');
    % Make legend
    ListOfLegendNames = {'MIGHT BE responsive genes'; 'KNOWN responsive genes'};
    legend(ListOfLegendNames);
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeast'); 
    hold off

    subplot(2,2,4)
    hist(ListOfListsToBePlot{8},20);
    h = findobj(gca,'Type','patch');
    hold on
    hist(ListOfListsToBePlot{6},20);
    % Color the 2 histograms differently
    set(findobj(gca,'Type','patch'),'FaceColor','r','EdgeColor','w');
    set(h,'FaceColor','b','EdgeColor','w');
    title('Treg - Cell 2'); xlabel('Test statistics Q (Ljung-Box test)'); ylabel('Number of genes per bin');
    % Make legend
    ListOfLegendNames = {'MIGHT BE responsive genes'; 'KNOWN responsive genes'};
    legend(ListOfLegendNames);
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeast'); 
    hold off
    %disp('arg')
    savefig('FILTERED_LjungBox_TestStatisticDistribution_SELECTED_Genes');
end
    
if MakeHistogramsOfLjungTestStatisticsALL_FILTERED == 1
    ListOfListsToBePlot = {};
    for TcellType = {'TeffALLgen', 'TregALLgen', 'Teff_MIGHT','Treg_MIGHT', 'Teff_KNOWN','Treg_KNOWN'}

        FileNameToBeLoaded = ['FILTERED_ExtractedCutData_' , TcellType{1}, '.mat'];
        load(FileNameToBeLoaded)
        for TableToBeUsed = {LjungBoxTestQandCL_Tcell1, LjungBoxTestQandCL_Tcell2}
            %disp(TableToBeUsed{1});
            %disp(TableToBeUsed{1}(:,1));
            ListOfListsToBePlot{end+1} = transpose(TableToBeUsed{1}(:,1));
        end

    end
    
    % Add a X2 distribtuion for 15 d.o.f. on top of it
    xForX2 = 0:0.1:120;
    NLags = 15;
    X2distribution = chi2pdf(xForX2,NLags);

    % Make Histograms of the distribution of the Test Statistics for all
    % Treg/Teff x cell1/cell2 x ALL genes
    figure

    subplot(2,2,1)
    histogram(ListOfListsToBePlot{1},150,'Normalization','pdf');
    hold on
    histogram(ListOfListsToBePlot{5},50,'Normalization','pdf');
    histogram(ListOfListsToBePlot{9},20,'Normalization','pdf');
    plot(xForX2,X2distribution)
    title('Teff - Cell 1'); xlabel('Test statistics Q (Ljung-Box test)'); ylabel('Probability density');
    % Make legend
    ListOfLegendNames = {'ALL FILTERED genes', 'FILTERED MIGHT BE responsive genes', 'FILTERED KNOWN responsive genes', '$$\chi^2$$ distribution, 15 d.o.f.'};
    legend(ListOfLegendNames,'Interpreter','latex');
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeast');
    hold off

    subplot(2,2,2)
    histogram(ListOfListsToBePlot{3},150,'Normalization','pdf');
    hold on
    histogram(ListOfListsToBePlot{7},50,'Normalization','pdf');
    histogram(ListOfListsToBePlot{11},20,'Normalization','pdf');
    plot(xForX2,X2distribution)
    title('Treg - Cell 1'); xlabel('Test statistics Q (Ljung-Box test)'); ylabel('Probability density');
    % Make legend    
    ListOfLegendNames = {'ALL FILTERED genes', 'FILTERED MIGHT BE responsive genes', 'FILTERED KNOWN responsive genes', '$$\chi^2$$ distribution, 15 d.o.f.'};
    legend(ListOfLegendNames,'Interpreter','latex');
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeast'); 
    hold off

    subplot(2,2,3)
    histogram(ListOfListsToBePlot{2},150,'Normalization','pdf');
    hold on
    histogram(ListOfListsToBePlot{6},50,'Normalization','pdf');
    histogram(ListOfListsToBePlot{10},20,'Normalization','pdf');
    plot(xForX2,X2distribution)
    title('Teff - Cell 2'); xlabel('Test statistics Q (Ljung-Box test)'); ylabel('Probability density');
    % Make legend
    ListOfLegendNames = {'ALL FILTERED genes', 'FILTERED MIGHT BE responsive genes', 'FILTERED KNOWN responsive genes', '$$\chi^2$$ distribution, 15 d.o.f.'};
    legend(ListOfLegendNames,'Interpreter','latex');
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeast'); 
    hold off

    subplot(2,2,4)
    histogram(ListOfListsToBePlot{4},150,'Normalization','pdf');
    hold on
    histogram(ListOfListsToBePlot{8},50,'Normalization','pdf');
    histogram(ListOfListsToBePlot{12},20,'Normalization','pdf');
    plot(xForX2,X2distribution)
    title('Treg - Cell 2'); xlabel('Test statistics Q (Ljung-Box test)'); ylabel('Probability density');
    % Make legend
    ListOfLegendNames = {'ALL FILTERED genes', 'FILTERED MIGHT BE responsive genes', 'FILTERED KNOWN responsive genes', '$$\chi^2$$ distribution, 15 d.o.f.'};
    legend(ListOfLegendNames,'Interpreter','latex');
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeast'); 
    hold off
    %disp('arg2')
    savefig('FILTERED_LjungBox_TestStatisticDistribution_ALL_Genes');
end