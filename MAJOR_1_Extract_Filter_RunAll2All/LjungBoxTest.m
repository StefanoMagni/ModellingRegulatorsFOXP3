
MakeIllustrativePlots = 1;
MakeHistogramsOfLjungTestStatisticsKNOWNandMIGHT = 1;
%ExtractGenesDataForComputingLjungTestStatistics = 0;
MakeHistogramsOfLjungTestStatisticsALL = 1;
LjungBoxTestOnPureWitheNoise = 1;

%%%%%%%%%% Make some plots to clarify meaning of quantities and tests used %%%%%%%%%%
if MakeIllustrativePlots == 1
    % Select a list and a couple of genes to be compared
    for MyGenesSelectionForNoiseFilteringTest = {{'ExtractedCutData_Teff_KNOWN.mat', 1, 13}, ...
            {'ExtractedCutData_Treg_MIGHT.mat',283,272}, {'ExtractedCutData_Treg_MIGHT.mat',208,55}}

        disp(MyGenesSelectionForNoiseFilteringTest)

        GenesListName = MyGenesSelectionForNoiseFilteringTest{1}{1};
        GoodGenePosition = MyGenesSelectionForNoiseFilteringTest{1}{2};
        BadGenePosition = MyGenesSelectionForNoiseFilteringTest{1}{3};

        load(GenesListName); % load('ExtractedCutData_Teff_MIGHT.mat');
        MyGoodGene = Tcell2DataTable(GoodGenePosition,:); % Tcell2DataTable(283,:); Tcell2DataTable(208,:);
        MyNoisyGene =Tcell2DataTable(BadGenePosition,:); % Tcell2DataTable(272,:); Tcell2DataTable(55,:);

        % Plot Gene Expression Rate for 1 noisy and one clearly responding gene
        figure
        plot(MyGoodGene)
        hold on
        plot(MyNoisyGene)
        hold off

        % Fast Fourier Transform (Jorge)
        figure
        plot(abs(fft(MyGoodGene-mean(MyGoodGene))));
        hold on
        plot(abs(fft(MyNoisyGene-mean(MyNoisyGene))));
        hold off

        % Ljung Box Test (Atte)
        figure
        subplot(2,1,1)
        autocorr(MyGoodGene)
        subplot(2,1,2)
        autocorr(MyNoisyGene)

        y = [];
        for i = 1:size(Tcell2DataTable,1) 
            MyGene = autocorr(Tcell2DataTable(i,:)); 
            x = 19*(19+2)*sum(MyGene(2:15).^2./(18:-1:5));
            y=[y,x];
        end
        disp(y);
    end
end

%%%%%%%%%% Now compute the test statistic Q for each gene of KNOWN/MIGHT %%%%%%%%%%
if MakeHistogramsOfLjungTestStatisticsKNOWNandMIGHT == 1
    NLags = 15;
    ListOfListsToBePlot = {};
    ListOfListsOfCL = {};
    for TcellType = {'Teff_KNOWN', 'Teff_MIGHT', 'Treg_KNOWN', 'Treg_MIGHT'}

        FileNameToBeLoaded = ['ExtractedCutData_' , TcellType{1}, '.mat'];
        load(FileNameToBeLoaded)

        for DataTableToBeUsed = {Tcell1DataTable, Tcell2DataTable}

            ListOfTestStatisticValues = [];
            %ListOfCLfromChiSquareDist = []
            for i = 1:size(DataTableToBeUsed{1},1) %18 
                DataPointsVector = DataTableToBeUsed{1}(i,:);
                TestStatistics = ComputeLjungBoxTestStatistics(DataPointsVector, NLags);
                ListOfTestStatisticValues=[ListOfTestStatisticValues,TestStatistics];
                %ConfidenceLevel = xxx
                %ListOfCLfromChiSquareDist=[ListOfCLfromChiSquareDist, ConfidenceLevel];
%                 MyGeneAutocorrelation = autocorr(DataTableToBeUsed{1}(i,:)); 
%                 TestStatistics = 19*(19+2)*sum(MyGeneAutocorrelation(2:15).^2./(18:-1:5));
%                 ListOfTestStatisticValues=[ListOfTestStatisticValues,TestStatistics];
            end

            % Use the cumulative distribution of the Chi square
            % distribution to compute the Confidence Level of considering a
            % gene's time-series as due to noise.
            NDegOfFreedom = NLags;
            ListOfConfidenceLevels = chi2cdf(ListOfTestStatisticValues,NDegOfFreedom);
            %disp(ListOfTestStatisticValues);
            %disp(ConfidenceLevel);
            
            ListOfTestStatisticValues
            [a,b] = min(ListOfTestStatisticValues)
            [c,d] = max(ListOfTestStatisticValues)

            ListOfListsToBePlot{end+1} = ListOfTestStatisticValues;
            ListOfListsOfCL{end+1} = ListOfConfidenceLevels;
        end
        % Save values of Ljung Box Test Statitstics Q into old datafiles
        LjungBoxTestQandCL_Tcell1 = transpose([ListOfListsToBePlot{end-1};ListOfListsOfCL{end-1}]);
        LjungBoxTestQandCL_Tcell2 = transpose([ListOfListsToBePlot{end};ListOfListsOfCL{end}]);
        save(FileNameToBeLoaded, 'LjungBoxTestQandCL_Tcell1','LjungBoxTestQandCL_Tcell2','-append');
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
    ListOfLegendNames = {'MIGHT BE responsive genes'; 'KNOWN responsive genes'};
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
    
    savefig('LjungBox_TestStatisticDistribution_SELECTED_Genes');
end


% %%%%%%%%%% Extract Data For Computing Ljung Test Statistics for ALL genes %%%%%%%%%%
% if ExtractGenesDataForComputingLjungTestStatistics == 1
%     % Teff
%     GenesDataExtractor(2, {}, 'Teff', 0)
%     % Treg
%     GenesDataExtractor(2, {}, 'Treg', 0)
% end


%%%%%%%%% Now compute the test statistics Q for ALL GENES!!! %%%%%%%%%
if MakeHistogramsOfLjungTestStatisticsALL == 1
    NLags = 15;
    ListOfListsToBePlot = {};
    ListOfListsOfCL = {};
    for TcellType = {'TeffALLgen', 'TregALLgen', 'Teff_MIGHT','Treg_MIGHT', 'Teff_KNOWN','Treg_KNOWN'}

        FileNameToBeLoaded = ['ExtractedCutData_' , TcellType{1}, '.mat'];
        load(FileNameToBeLoaded)

        for DataTableToBeUsed = {Tcell1DataTable, Tcell2DataTable}

            ListOfTestStatisticValues = [];
            for i = 1:size(DataTableToBeUsed{1},1)%18 
                DataPointsVector = DataTableToBeUsed{1}(i,:);
                TestStatistics = ComputeLjungBoxTestStatistics(DataPointsVector, NLags);
                ListOfTestStatisticValues=[ListOfTestStatisticValues,TestStatistics];
            end
            
            % Use the cumulative distribution of the Chi square
            % distribution to compute the Confidence Level of considering a
            % gene's time-series as due to noise.
            NDegOfFreedom = NLags;
            ListOfConfidenceLevels = chi2cdf(ListOfTestStatisticValues,NDegOfFreedom);
            %disp(ListOfTestStatisticValues);
            %disp(ConfidenceLevel);
            
            [a,b] = min(ListOfTestStatisticValues)
            [c,d] = max(ListOfTestStatisticValues)

            ListOfListsToBePlot{end+1} = ListOfTestStatisticValues;
            ListOfListsOfCL{end+1} = ListOfConfidenceLevels;
        end
        % Save values of Ljung Box Test Statitstics Q into old datafiles
        LjungBoxTestQandCL_Tcell1 = transpose([ListOfListsToBePlot{end-1};ListOfListsOfCL{end-1}]);
        LjungBoxTestQandCL_Tcell2 = transpose([ListOfListsToBePlot{end};ListOfListsOfCL{end}]);
        save(FileNameToBeLoaded, 'LjungBoxTestQandCL_Tcell1','LjungBoxTestQandCL_Tcell2','-append');
    end

    % Add a X2 distribtuion for 15 d.o.f. on top of it
    xForX2 = 0:0.1:120;
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
    ListOfLegendNames = {'ALL genes', 'MIGHT BE responsive genes', 'KNOWN responsive genes', '$$\chi^2$$ distribution, 15 d.o.f.'};
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
    ListOfLegendNames = {'ALL genes', 'MIGHT BE responsive genes', 'KNOWN responsive genes', '$$\chi^2$$ distribution, 15 d.o.f.'};
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
    ListOfLegendNames = {'ALL genes', 'MIGHT BE responsive genes', 'KNOWN responsive genes', '$$\chi^2$$ distribution, 15 d.o.f.'};
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
    ListOfLegendNames = {'ALL genes', 'MIGHT BE responsive genes', 'KNOWN responsive genes', '$$\chi^2$$ distribution, 15 d.o.f.'};
    legend(ListOfLegendNames,'Interpreter','latex');
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeast'); 
    hold off
    
    savefig('LjungBox_TestStatisticDistribution_ALL_Genes');
end

if LjungBoxTestOnPureWitheNoise == 1
    figure
    hold on
    for NLags = [5,10,18]
        % Create fake data with white noise only
        mu = 0;
        sigma = 500;
        ListOfTestStatisticsOnWhiteNoise = [];
        for i = 1:1:10000
            MyWhiteNoiseTimeSeriesData = normrnd(mu,sigma,[1,19]);
            % And compute the test statistics for Ljung Box test
            Q = ComputeLjungBoxTestStatistics(MyWhiteNoiseTimeSeriesData,NLags);
            ListOfTestStatisticsOnWhiteNoise = [ListOfTestStatisticsOnWhiteNoise, Q];
        end
        % Add a X2 distribtuion for NLags d.o.f. on top of it
        xForX2 = 0:0.1:100;
        X2distribution = chi2pdf(xForX2,NLags);
        % Make Histograms of the distribution of the Test Statistics for all
        histogram(ListOfTestStatisticsOnWhiteNoise,100,'Normalization','pdf');
        plot(xForX2,X2distribution)
    end
    title('Ljung Box Test Statistic for White Noise time-series'); 
    xlabel('Test statistics Q (Ljung-Box test)'); 
    ylabel('Probability density');
    xlim([0 50]);
    % Make legend
    ListOfLegendNames = {'10000 white noise samples, NLag = 5', '$$\chi^2$$ distribution, 5 d.o.f.', ...
        '10000 white noise samples, NLag = 10', '$$\chi^2$$ distribution, 10 d.o.f.', ...
        '10000 white noise samples, NLag = 18', '$$\chi^2$$ distribution, 18 d.o.f.'};
    legend(ListOfLegendNames,'Interpreter','latex');
    lh=findall(gcf,'tag','legend');
    set(lh,'location','northeast'); 
    hold off
    savefig('LjungBox_TestStatisticDistribution_WHITE_NOISE');
end

%%%%%%%%%% DEFINE FUNCTION TO COMPUTE LJUNG BOX TEST STATISTIC %%%%%%%%%%
function [TestStatistics] = ComputeLjungBoxTestStatistics(DataPointsVector, NLags)
Ndatapoins = numel(DataPointsVector);
MyGeneAutocorrelation = autocorr(DataPointsVector); 
TestStatistics = Ndatapoins * (Ndatapoins+2) * ...
    sum(MyGeneAutocorrelation(2:NLags+1).^2 ./ ...
    ((Ndatapoins-1):-1:(Ndatapoins-NLags)));
end

