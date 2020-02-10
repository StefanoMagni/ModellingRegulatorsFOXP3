% Author: Stefano Magni, created on 26/12/2016, contact: stefano.magni@uni.lu
% Script to extract desired portions of datafiles (gene expression data)
% to reproduce plots in He 2012 paper figure 2

%clear, close all, clc;

ListOfGeneNames = {'IL2';'IL4';'CCL20';'IL5';'CSF2';'FOXP3';'IL13';'IFNG';'LRRC32'};

GenesDataExtractor(1, ListOfGeneNames, 'Treg', 1);
GenesDataExtractor(1, ListOfGeneNames, 'Teff', 1);

figure;

for TcellTYPEloop = {'Treg','Teff'}

    CurrentFileName = ['ExtractedCutData_', TcellTYPEloop{1}, 'Custom.mat'];
    load(CurrentFileName);
    
    if TcellTYPEloop{1} == 'Treg'
        ColorList = {'o-blue';'s-magenta'};
    elseif TcellTYPEloop{1} == 'Teff'
        ColorList = {'^-Yellow';'x-cyan'};
    end

    % PLOT THE DATAPOINTS FOR ALL GENES AS FUNC OF t
    TimePointsList = [0.]; % (minutes)
    for i = 1:(size(Tcell1DataTable,2)-1)
        TimePointsList = [ TimePointsList, TimePointsList(i) + 20.];
    end
    
    j = 1;
    for k = 1:9
        % Make subplot k
        Dj = cell2mat(ListOfGeneMatches(k,2));
        subplot(3,3,k);
        for jj = j:j+Dj-1
            plot(TimePointsList,Tcell1DataTable(jj,:),ColorList{1}); hold on;
            plot(TimePointsList,Tcell2DataTable(jj,:),ColorList{2}); hold on;
        end
        j = j+Dj;
        title(ListOfGeneNames(k)); xlabel('time (min)'); ylabel('mRNA concentration (a.u.?)');
        hold on;
    end
    hold on;
end

% Make legend
ListOfLegendNames = ['Treg1';'Treg2';'Teff1';'Teff2'];
legend(ListOfLegendNames);
lh=findall(gcf,'tag','legend');
set(lh,'location','northeastoutside'); 
hold off;

% Make general plot title
annotation('textbox', [0 0.9 1 0.1], ...
'String', ['Gene expression time-series as in Fig. 2 of He 2012'], ...
'EdgeColor', 'none', ...
'HorizontalAlignment', 'center');

% Give a proper name to the figure and save it
FigName = sprintf(['GenesExpression_ReproducingFig2He2012']);
savefig(FigName);

