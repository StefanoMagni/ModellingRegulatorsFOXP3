%%%% Code to run PCA (SVD) on Treg data
clear; close all; clc;

Do_PCA_SVD = 0;
NetworkInference = 0;
PlotBiograph = 1;
PileUpProbes = 0; % PROBABLY WRONG !!! DO NOT USE THIS !!!
% DANGER!!! THIS STILL NEED TO BE FIXED, FOR NOW DCGAIN IS 
% TAKEN RANDOMLY, WHILE THE DC GAIN WHEN PILING UP PROBES SHOULD BE TAKEN AS 
% THE ONE CORRESPONDING TO HIGHEST FITNESS LINK, AND IT SHOULD BE TESTED IF 
% DIFFERENT BROBES PROPVIDES DIFFERENT ACTIVATION/REPRESSION, AND IDEALLY 
% IT SHOULD BE INDICATED THAT A NODE COMES INDEED FROM MORE THAN 1 PROBE!!!

if Do_PCA_SVD == 1;
    % M = [1 1 0; 2 2 0; (1/1.414) (1/1.414) 0; 1 0 0]
    % 
    % [U,S,V] = svd(M)
    % 
    % scatter(M(:,1),M(:,2));
    % 
    % % V2 = pca(M) % This could be run alternatively to get only the PCs
    % 
    % mapcaplot(M)

    %%%%%%%%%%%%%%%%%%%%

    % load('FILTERED_ExtractedCutData_Teff_MIGHT.mat', 'Tcell1DataTable', 'Tcell2DataTable')
    % DataTeff1 = Tcell1DataTable;
    % DataTeff2 = Tcell2DataTable;
    % 
    % load('FILTERED_ExtractedCutData_Treg_MIGHT.mat', 'Tcell1DataTable', 'Tcell2DataTable')
    % DataTreg1 = Tcell1DataTable;
    % DataTreg2 = Tcell2DataTable;

    % mapcaplot(DataTeff2'); 
    % The following alternative can be used as well:
    % [coeff,score] = pca(DataTeff2');
    % scatter(score(:,1),score(:,2));

    load('FILTERED_ExtractedCutData_TeffALLgen.mat', 'Tcell1DataTable', 'Tcell2DataTable')
    DataTeff1ALL = Tcell1DataTable;
    DataTeff2ALL = Tcell2DataTable;

    load('FILTERED_ExtractedCutData_TregALLgen.mat', 'Tcell1DataTable', 'Tcell2DataTable')
    DataTreg1ALL = Tcell1DataTable;
    DataTreg2ALL = Tcell2DataTable;

    % mapcaplot(DataTreg2ALL');
    MyColors = 20*[0:18];
    pointsize = 100;

    [coeffTeff1,scoreTeff1,latentTeff1,tsquaredTeff1,explainedTeff1,muTeff1] = pca(DataTeff1ALL');
    [coeffTeff2,scoreTeff2,latentTeff2,tsquaredTeff2,explainedTeff2,muTeff2] = pca(DataTeff2ALL');
    [coeffTreg1,scoreTreg1,latentTreg1,tsquaredTreg1,explainedTreg1,muTreg1] = pca(DataTreg1ALL');
    [coeffTreg2,scoreTreg2,latentTreg2,tsquaredTreg2,explainedTreg2,muTreg2] = pca(DataTreg2ALL');

    %%%%%%%%%
    figure;

    subplot(2,2,1)
    explainedTeff1MOD = [explainedTeff1(1); explainedTeff1(2); explainedTeff1(3); sum(explainedTeff1(4:18))];
    pie(explainedTeff1MOD);
    legend(['PC1';'PC2';'PC3';'...'], 'Location','NorthWest');
    title('Teff 1');

    subplot(2,2,3)
    explainedTeff2MOD = [explainedTeff2(1); explainedTeff2(2); explainedTeff2(3); sum(explainedTeff2(4:18))];
    pie(explainedTeff2MOD);
    title('Teff 2');

    subplot(2,2,2)
    explainedTreg1MOD = [explainedTreg1(1); explainedTreg1(2); explainedTreg1(3); sum(explainedTreg1(4:18))];
    pie(explainedTreg1MOD);
    title('Treg 1');

    subplot(2,2,4)
    explainedTreg2MOD = [explainedTreg2(1); explainedTreg2(2); explainedTreg2(3); sum(explainedTreg2(4:18))];
    pie(explainedTreg2MOD);
    title('Treg 2');

    %%%%%%%%%
    figure;

    subplot(2,2,1)
    scatter(scoreTeff1(:,1),scoreTeff1(:,2), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Teff Cell 1');
    xlabel('PC 1');
    ylabel('PC 2');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,3)
    scatter(scoreTeff2(:,1),scoreTeff2(:,2), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Teff Cell 2');
    xlabel('PC 1');
    ylabel('PC 2');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,2)
    scatter(scoreTreg1(:,1),scoreTreg1(:,2), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Treg Cell 1');
    xlabel('PC 1');
    ylabel('PC 2');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,4)
    scatter(scoreTreg2(:,1),scoreTreg2(:,2), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Treg Cell 2');
    xlabel('PC 1');
    ylabel('PC 2');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;

    subplot(2,2,1)
    scatter(scoreTeff1(:,1),scoreTeff1(:,3), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Teff Cell 1');
    xlabel('PC 1');
    ylabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,3)
    scatter(scoreTeff2(:,1),scoreTeff2(:,3), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Teff Cell 2');
    xlabel('PC 1');
    ylabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,2)
    scatter(scoreTreg1(:,1),scoreTreg1(:,3), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Treg Cell 1');
    xlabel('PC 1');
    ylabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,4)
    scatter(scoreTreg2(:,1),scoreTreg2(:,3), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Treg Cell 2');
    xlabel('PC 1');
    ylabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;

    subplot(2,2,1)
    scatter3(scoreTeff1(:,1),scoreTeff1(:,2),scoreTeff1(:,3), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Teff Cell 1');
    xlabel('PC 1');
    ylabel('PC 2');
    zlabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,3)
    scatter3(scoreTeff2(:,1),scoreTeff2(:,2),scoreTeff2(:,3), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Teff Cell 2');
    xlabel('PC 1');
    ylabel('PC 2');
    zlabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,2)
    scatter3(scoreTreg1(:,1),scoreTreg1(:,2),scoreTreg1(:,3), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Treg Cell 1');
    xlabel('PC 1');
    ylabel('PC 2');
    zlabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,4)
    scatter3(scoreTreg2(:,1),scoreTreg2(:,2),scoreTreg2(:,3), pointsize, MyColors, 'MarkerFaceColor', 'flat');
    title('Treg Cell 2');
    xlabel('PC 1');
    ylabel('PC 2');
    zlabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;

    subplot(2,2,1)
    plot3(scoreTeff1(:,1),scoreTeff1(:,2),scoreTeff1(:,3), 'Color', 'b', 'LineWidth', 1, 'Marker', 'o');
    title('Teff Cell 1');
    xlabel('PC 1');
    ylabel('PC 2');
    zlabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,3)
    plot3(scoreTeff2(:,1),scoreTeff2(:,2),scoreTeff2(:,3), 'Color', 'r', 'LineWidth', 1, 'Marker', 'o');
    title('Teff Cell 2');
    xlabel('PC 1');
    ylabel('PC 2');
    zlabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,2)
    plot3(scoreTreg1(:,1),scoreTreg1(:,2),scoreTreg1(:,3), 'Color', 'g', 'LineWidth', 1, 'Marker', 'o');
    title('Treg Cell 1');
    xlabel('PC 1');
    ylabel('PC 2');
    zlabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';

    subplot(2,2,4)
    plot3(scoreTreg2(:,1),scoreTreg2(:,2),scoreTreg2(:,3), 'Color', 'black', 'LineWidth', 1, 'Marker', 'o');
    title('Treg Cell 2');
    xlabel('PC 1');
    ylabel('PC 2');
    zlabel('PC 3');
    c = colorbar;
    c.Label.String = 'Time (minutes)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Ue1,Se1,Ve1] = svd(DataTeff1ALL);
    [Ur1,Sr1,Vr1] = svd(DataTreg1ALL);
    [Ue2,Se2,Ve2] = svd(DataTeff2ALL);
    [Ur2,Sr2,Vr2] = svd(DataTreg2ALL);

    % figure;
    % for i = [1:3]
    %     plot(20*[0:18],Ve2(:,i)); hold on;
    % end
    % title('PC1 PC2 PC3 basis vectors - Teff 2');
    % xlabel('time (minutes)');
    % ylabel('PC');
    % 
    % %%%%%%%%%%%%%%%%%%%%
    % 
    % figure;
    % for i = [1:3]
    %     plot(20*[0:18],Vr2(:,i)); hold on;
    % end
    % title('PC1 PC2 PC3 basis vectors - Treg 2');
    % xlabel('time (minutes)');
    % ylabel('PC');

    %%%%%%%%%%%%%%%%%%%%%%%%%%% PCs BASIS VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;

    subplot(2,2,1)
    for i = [1:3]
        plot(20*[0:18],Ve1(:,i),'LineWidth',1.6); hold on;
    end
    plot(20*[0:18],zeros(19,1), 'Color','black');
    title('Teff 1');
    xlim([0 20*18])
    ylim([-0.5 1])
    xlabel('time (minutes)');
    ylabel('PC');
    legend(['PC1'; 'PC2'; 'PC3'],'Location','NorthWest');

    subplot(2,2,2)
    for i = [1:3]
        plot(20*[0:18],Vr1(:,i),'LineWidth',1.6); hold on;
    end
    plot(20*[0:18],zeros(19,1), 'Color','black');
    title('Treg 1');
    xlim([0 20*18])
    ylim([-0.5 1])
    xlabel('time (minutes)');
    ylabel('PC');
    legend(['PC1'; 'PC2'; 'PC3'],'Location','NorthWest');

    subplot(2,2,3)
    for i = [1:3]
        plot(20*[0:18],Ve2(:,i),'LineWidth',1.6); hold on;
    end
    plot(20*[0:18],zeros(19,1), 'Color','black');
    title('Teff 2');
    xlim([0 20*18])
    ylim([-0.5 1])
    xlabel('time (minutes)');
    ylabel('PC');
    legend(['PC1'; 'PC2'; 'PC3'],'Location','NorthWest');

    subplot(2,2,4)
    for i = [1:3]
        plot(20*[0:18],Vr2(:,i),'LineWidth',1.6); hold on;
    end
    plot(20*[0:18],zeros(19,1), 'Color','black');
    title('Treg 2');
    xlim([0 20*18])
    ylim([-0.5 1])
    xlabel('time (minutes)');
    ylabel('PC');
    legend(['PC1'; 'PC2'; 'PC3'],'Location','NorthWest');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load('FOXP3_TregCell12_Probes23.mat', 'FOXP3_Probe2_Treg12')
    load('FOXP3_TregCell12_Probes23.mat', 'FOXP3_Probe3_Treg12')
    CumultaiveFOXP3exprTreg1 = (FOXP3_Probe2_Treg12(1,:) + FOXP3_Probe3_Treg12(1,:)) / 2;
    CumultaiveFOXP3exprTreg2 = (FOXP3_Probe2_Treg12(2,:) + FOXP3_Probe3_Treg12(2,:)) / 2;
    figure;
    plot(20*[0:18],CumultaiveFOXP3exprTreg1); hold on;
    plot(20*[0:18],CumultaiveFOXP3exprTreg2); hold off;

    figure;

    subplot(1,2,1);
    scatter(scoreTreg1(:,1),scoreTreg1(:,2), pointsize, CumultaiveFOXP3exprTreg1, 'MarkerFaceColor', 'flat');
    title('Treg Cell 1');
    xlabel('PC 1');
    ylabel('PC 2');
    c = colorbar;
    c.Label.String = 'FOXP3 Epression';
    a = [0:18]'; b = num2str(a); c = cellstr(b);
    dx = 1000; dy = 1000; % displacement so the text does not overlay the data points
    text(scoreTreg1(:,1)+dx, scoreTreg1(:,2)+dy, c);

    subplot(1,2,2);
    scatter(scoreTreg2(:,1),scoreTreg2(:,2), pointsize, CumultaiveFOXP3exprTreg2, 'MarkerFaceColor', 'flat');
    title('Treg Cell 2');
    xlabel('PC 1');
    ylabel('PC 2');
    c = colorbar;
    c.Label.String = 'FOXP3 Epression';
    a = [0:18]'; b = num2str(a); c = cellstr(b);
    dx = 1000; dy = 1000; % displacement so the text does not overlay the data points
    text(scoreTreg2(:,1)+dx, scoreTreg2(:,2)+dy, c);
    %%%%%%%%%%%%%%%%%%%%

    load('IKZF4_TregCell12_Probes234.mat', 'IKZF4_Probe2_Treg12')
    load('IKZF4_TregCell12_Probes234.mat', 'IKZF4_Probe3_Treg12')
    load('IKZF4_TregCell12_Probes234.mat', 'IKZF4_Probe4_Treg12')
    CumultaiveIKZF4exprTreg1 = (IKZF4_Probe2_Treg12(1,:) + IKZF4_Probe3_Treg12(1,:) + IKZF4_Probe4_Treg12(1,:)) / 3;
    CumultaiveIKZF4exprTreg2 = (IKZF4_Probe2_Treg12(2,:) + IKZF4_Probe3_Treg12(2,:) + IKZF4_Probe4_Treg12(2,:)) / 3;
    figure;
    plot(20*[0:18],CumultaiveIKZF4exprTreg1); hold on;
    plot(20*[0:18],CumultaiveIKZF4exprTreg2); hold off;

    figure;

    subplot(1,2,1);
    scatter(scoreTreg1(:,1),scoreTreg1(:,2), pointsize, CumultaiveIKZF4exprTreg1, 'MarkerFaceColor', 'flat');
    title('Treg Cell 1');
    xlabel('PC 1');
    ylabel('PC 2');
    c = colorbar;
    c.Label.String = 'IKZF4 Epression';
    a = [0:18]'; b = num2str(a); c = cellstr(b);
    dx = 1000; dy = 1000; % displacement so the text does not overlay the data points
    text(scoreTreg1(:,1)+dx, scoreTreg1(:,2)+dy, c);

    subplot(1,2,2);
    scatter(scoreTreg2(:,1),scoreTreg2(:,2), pointsize, CumultaiveIKZF4exprTreg2, 'MarkerFaceColor', 'flat');
    title('Treg Cell 2');
    xlabel('PC 1');
    ylabel('PC 2');
    c = colorbar;
    c.Label.String = 'IKZF4 Epression';
    a = [0:18]'; b = num2str(a); c = cellstr(b);
    dx = 1000; dy = 1000; % displacement so the text does not overlay the data points
    text(scoreTreg2(:,1)+dx, scoreTreg2(:,2)+dy, c);

    %%%%%%%%%%%%%%%%%%%%

    load('ExtractedCutData_TeffALLgen.mat', 'Tcell1DataTable');
    load('ExtractedCutData_TeffALLgen.mat', 'Tcell2DataTable');
    DataTeff1ALLunfiltered = Tcell1DataTable;
    DataTeff2ALLunfiltered = Tcell2DataTable;
    load('ExtractedCutData_TregALLgen.mat', 'Tcell1DataTable');
    load('ExtractedCutData_TregALLgen.mat', 'Tcell2DataTable');
    DataTreg1ALLunfiltered = Tcell1DataTable;
    DataTreg2ALLunfiltered = Tcell2DataTable;

    AllUnfilteredData = [DataTeff1ALLunfiltered, DataTeff2ALLunfiltered, ...
                         DataTreg1ALLunfiltered, DataTreg2ALLunfiltered];

    MyNewColors = [ones(1,19), 2*ones(1,19), 3*ones(1,19), 4*ones(1,19)];
    [coeff,score] = pca(AllUnfilteredData');

    figure;
    scatter(score(:,1),score(:,2), pointsize, MyNewColors, 'MarkerFaceColor', 'flat');
    title('Teff 1, Teff 2, Treg 1, Treg 2');
    xlabel('PC 1');
    ylabel('PC 2');
    colorbar('Ticks',[1,2,3,4],...
             'TickLabels',{'Teff 1','Teff 2','Treg 1','Treg 2'})

    figure;
    scatter3(score(:,1),score(:,2), score(:,3), pointsize, MyNewColors, 'MarkerFaceColor', 'flat');
    title('Teff 1, Teff 2, Treg 1, Treg 2');
    xlabel('PC 1');
    ylabel('PC 2');
    zlabel('PC 3');
    colorbar('Ticks',[1,2,3,4],...
             'TickLabels',{'Teff 1','Teff 2','Treg 1','Treg 2'})

    figure;
    plot3(score(1:19,1),score(1:19,2), score(1:19,3), 'Color', 'b', 'LineWidth', 1, 'Marker', 'o'); hold on;
    plot3(score(20:38,1),score(20:38,2), score(20:38,3), 'Color', 'r', 'LineWidth', 1, 'Marker', 'o'); hold on;
    plot3(score(39:57,1),score(39:57,2), score(39:57,3), 'Color', 'c', 'LineWidth', 1, 'Marker', 'o'); hold on;
    plot3(score(58:76,1),score(58:76,2), score(58:76,3), 'Color', 'g', 'LineWidth', 1, 'Marker', 'o'); hold on;
    legend(['Teff 1'; 'Teff 2'; 'Treg 1'; 'Treg 2']);
    xlabel('PC 1');
    ylabel('PC 2');
    zlabel('PC 3');

    % %%%%%% COMPUTE AVERAGE GENE GENE CORRELATION %%%%%%%
    % GeneGeneCorrelationTreg2 = corr(DataTreg2ALL(:,1)); % Pearson gene-gene correlation matrix
    % MeanGeneGeneCorrelationtreg2 = mean(mean(GeneGeneCorrelationTreg2));

    %%%%% COMPUTE EUCLIDIAN DISTANCE BETWEEN EACH 2 COUPLES OF VECTORS %%%%%
    TimeSeriesEuclidianDistancesTeff1 = [];
    for i = 1: 18
        EuclidianDistance = norm(DataTeff1ALL(:,i)-DataTeff1ALL(:,i+1));
        TimeSeriesEuclidianDistancesTeff1 = [TimeSeriesEuclidianDistancesTeff1, ...
                                             EuclidianDistance];
    end

    TimeSeriesEuclidianDistancesTeff2 = [];
    for i = 1: 18
        EuclidianDistance = norm(DataTeff2ALL(:,i)-DataTeff2ALL(:,i+1));
        TimeSeriesEuclidianDistancesTeff2 = [TimeSeriesEuclidianDistancesTeff2, ...
                                             EuclidianDistance];
    end

    TimeSeriesEuclidianDistancesTreg1 = [];
    for i = 1: 18
        EuclidianDistance = norm(DataTreg1ALL(:,i)-DataTreg1ALL(:,i+1));
        TimeSeriesEuclidianDistancesTreg1 = [TimeSeriesEuclidianDistancesTreg1, ...
                                             EuclidianDistance];
    end

    TimeSeriesEuclidianDistancesTreg2 = [];
    for i = 1: 18
        EuclidianDistance = norm(DataTreg2ALL(:,i)-DataTreg2ALL(:,i+1));
        TimeSeriesEuclidianDistancesTreg2 = [TimeSeriesEuclidianDistancesTreg2, ...
                                             EuclidianDistance];
    end

    % figure;
    % subplot(2,2,1); plot(TimeSeriesEuclidianDistancesTeff1); title('Teff 1');
    % subplot(2,2,2); plot(TimeSeriesEuclidianDistancesTreg1); title('Treg 1');
    % subplot(2,2,3); plot(TimeSeriesEuclidianDistancesTeff2); title('Teff 2');
    % subplot(2,2,4); plot(TimeSeriesEuclidianDistancesTreg2); title('Treg 2');

    figure;
    subplot(1,2,1); 
    plot(TimeSeriesEuclidianDistancesTeff1); hold on; 
    plot(TimeSeriesEuclidianDistancesTeff2); hold off;
    xlabel('Couple of consecutive timepoints');
    ylabel('Euclidian Distance');
    legend(['Teff 1'; 'Teff 2']);
    subplot(1,2,2); 
    plot(TimeSeriesEuclidianDistancesTreg1); hold on;
    plot(TimeSeriesEuclidianDistancesTreg2); hold off;
    xlabel('Couple of consecutive timepoints');
    ylabel('Euclidian Distance');
    legend(['Treg 1'; 'Treg 2']);

    %%%%%% COMPUTE CORRELATION BETWEEN NEIGHBOURS TIME VECTORS %%%%%

    TimeTimeCorrelationTeff1 = corr(DataTeff1ALL); % Pearson gene-gene correlation matrix
    TimeTimeNeighboursCorrelationTeff1 = [];
    for i = 1:18
        TimeTimeNeighboursCorrelationTeff1 = [TimeTimeNeighboursCorrelationTeff1, ...
            TimeTimeCorrelationTeff1(i,i+1)];
    end

    TimeTimeCorrelationTeff2 = corr(DataTeff2ALL); % Pearson gene-gene correlation matrix
    TimeTimeNeighboursCorrelationTeff2 = [];
    for i = 1:18
        TimeTimeNeighboursCorrelationTeff2 = [TimeTimeNeighboursCorrelationTeff2, ...
            TimeTimeCorrelationTeff2(i,i+1)];
    end

    TimeTimeCorrelationTreg1 = corr(DataTreg1ALL); % Pearson gene-gene correlation matrix
    TimeTimeNeighboursCorrelationTreg1 = [];
    for i = 1:18
        TimeTimeNeighboursCorrelationTreg1 = [TimeTimeNeighboursCorrelationTreg1, ...
            TimeTimeCorrelationTreg1(i,i+1)];
    end

    TimeTimeCorrelationTreg2 = corr(DataTreg2ALL); % Pearson gene-gene correlation matrix
    TimeTimeNeighboursCorrelationTreg2 = [];
    for i = 1:18
        TimeTimeNeighboursCorrelationTreg2 = [TimeTimeNeighboursCorrelationTreg2,
            TimeTimeCorrelationTreg2(i,i+1)];
    end

    figure;
    title('Time-Time correlation')
    subplot(1,2,1); 
    plot(TimeTimeNeighboursCorrelationTeff1); hold on; 
    plot(TimeTimeNeighboursCorrelationTeff2); hold off;
    xlabel('Couple of consecutive timepoints');
    ylabel('Gene-Gene Pearson Correlation');
    legend(['Teff 1'; 'Teff 2']);
    subplot(1,2,2); 
    plot(TimeTimeNeighboursCorrelationTreg1); hold on;
    plot(TimeTimeNeighboursCorrelationTreg2); hold off;
    xlabel('Couple of consecutive timepoints');
    ylabel('Gene-Gene Pearson Correlation');
    legend(['Treg 1'; 'Treg 2']);
end

% OBSERVATIONS/CONCLUSIONS:
% 1) Teff shows nicely continuous time evo
% 2) Treg insteead shows 3 clear phases, corresponding in Treg 2 to ...
%    FOXP3 highly expressed
% 3) FOXP3 higly expressed brings back to the region top right
% 4) Also IKZF4 correspond, even more nicely, to the split in 3 regions
% 5) The basis functions of Treg2 shows nicely 3 phases, unlike Teff
% 6) The PC1 base funct seems to be common among Treg and Teff, should it be due to
% the processes of Immune response which are common between the 2 genes?
% 7) 6h --> not yet back to normal!
% 8) 3-phase picture clear in Treg 2, less clear in Treg 1, absent in Teff
% 1/2, but synchro with FOXP3  only in Treg2, not in Treg1!!!
% 9) FOXP# seems to be bringing to phase II! But it seems that is not the
% only/main, as when cells enter phase 2 this gene is indeed not yet highly
% expressed. IKZF4?
% 10) If PC1 is the same between all 4 cells, and if it is a constant,
% may be biologically it represents the BASAL GENE EXPRESSION which would
% have been there without the stimulation mimiking the virus. 
% These genes are all expressed (or not expressed) with the same intensity 
% across different timepoints.
% Then movements of time-points along the direction of PC1 could simply
% represent that this basal methabolism is relatively more or less
% important (relative gene expression) with respect to the response to the
% stimultation.

%%%%%%%%%%%%%%%%%%%% Network inference (from Laurent's script) %%%%%%%%%%%%%%%%%%%%
if NetworkInference == 1
    
    load('FILTERED_ExtractedCutData_Treg_KNOWN.mat')

    % All2all Network Inference Options
    method = 'ATA';
    order = 1; %Order of system for ATA
    init_opt = 'n4sid'; %Choice between {'iv','svf','gpmf','n4sid','all'}. See doc "tfestoptions"

    opt = tfestOptions('InitMethod',init_opt); %Option for system ID ATA initial estimate technique. 

    for phase = {'__I','_II','III','ALL'}
        for CellIndex = 1:2
            %Run network inference technique (Do not consider self regulation)

            if CellIndex == 1
                data = Tcell1DataTable;
                %data = sscanf(sprintf('%s*', Tcell1DataTable), '%f*');
            elseif CellIndex == 2
                data = Tcell2DataTable;
            end

            if phase{1} == 'ALL'
                dataCUT = data;
            elseif phase{1} == '__I'
                dataCUT = data(:,1:13);
            elseif phase{1} == '_II'
                dataCUT = data(:,4:16);
            elseif phase{1} == 'III'
                dataCUT = data(:,7:19);
            end

            [fitness, models, IOs, AICs,dcgains] = just_tfest(order, 1, dataCUT, opt); %Only mRNA
            if CellIndex == 1
                mRNA_names = TableOfExtractedGeneNameAndIDs_cell1(:,4);
            elseif CellIndex == 2
                mRNA_names = TableOfExtractedGeneNameAndIDs_cell2(:,4);
            end
            ranks = ranking_list(fitness,mRNA_names,AICs);
            FileName = sprintf(['3STEPS_',phase{1},'_Results_Treg_KNOWN_Cell' num2str(CellIndex) '.mat']);
            save(FileName, 'fitness', 'models', 'IOs', 'AICs', 'ranks', 'mRNA_names', 'dcgains');
        end
    end
end
%%%%%%%%%%%%%%%%%%%% Plotting Biograph of the reconstructed networks %%%%%%%%%%%%%%%%%%%
    
if PlotBiograph == 1
    %load('FILTERED_ExtractedCutData_Treg_KNOWN.mat')
    load('FILTERED_ExtractedCutData_Treg_KNOWN.mat', 'TableOfExtractedGeneNameAndIDs_cell1');
    load('FILTERED_ExtractedCutData_Treg_KNOWN.mat', 'TableOfExtractedGeneNameAndIDs_cell2');
    
    % Network visual representation using biograph
    threshold = 70; % Treshold for All to All fitness score

    for phase = {'__I','_II','III','ALL'} 
        for CellIndex = 1:2%2
%             if PileUpProbes == 1
%                 if CellIndex == 1                    
%                     ResultsFileNameToBeLoaded = ['PiledUpProbes_3STEPS_',phase{1},'_Results_Treg_KNOWN_Cell1.mat'];
%                 elseif CellIndex == 2
%                     ResultsFileNameToBeLoaded = ['PiledUpProbes_3STEPS_',phase{1},'_Results_Treg_KNOWN_Cell2.mat'];
%                 end
%             end
%             if PileUpProbes == 0
            if CellIndex == 1                    
                ResultsFileNameToBeLoaded = ['3STEPS_',phase{1},'_Results_Treg_KNOWN_Cell1.mat'];
            elseif CellIndex == 2
                ResultsFileNameToBeLoaded = ['3STEPS_',phase{1},'_Results_Treg_KNOWN_Cell2.mat'];
            end
%             end
            load(ResultsFileNameToBeLoaded);
            
            % Chose table according to filtering / no filtering
            if CellIndex == 1
                MyTable = TableOfExtractedGeneNameAndIDs_cell1;
            elseif CellIndex == 2
                MyTable = TableOfExtractedGeneNameAndIDs_cell2;
            end

            if not(PileUpProbes)
                % Format Genes Names for Plotting
                mRNAs_names = {MyTable{1,4}};
                for i = 2:1:numel(MyTable(:,4))
                   k = 1;
                   HasTwins = 0;
                   TentativeNodeName = MyTable{i,4};
                    for j = 1:(i-1)
                        if strcmp(TentativeNodeName,mRNA_names{j})
                            k = k + 1;
                            HasTwins = 1;
                        end
                    end
                    if HasTwins
                        TentativeNodeName = [sprintf(TentativeNodeName), '-', num2str(k)];
                    else
                        TentativeNodeName = [sprintf(TentativeNodeName)];
                    end
                    mRNAs_names{i} = TentativeNodeName;
                    
                    % Make a Bio Graph representing the reconstracted gene regulatory network
                end
                fprintf('Plotting graph of reconstructed network...\n');
                FigName = sprintf(['3STEPS_',phase{1},'_GeneRegulatoryNetwork_Treg_KNOWN_Cell', num2str(CellIndex), '_Treg_FILTERED']);
                plot_biograph(fitness,threshold,mRNAs_names,FigName,dcgains);
                
            elseif PileUpProbes
                [CUTmRNA_names, CutCutNewMODfitness, CutCututdcgains] = PileUpProbesToOneNode( ...
                    MyTable, mRNA_names, fitness, dcgains)
                fprintf('Plotting graph of reconstructed network...\n');
                FigName = sprintf(['3STEPS_PiledUpProbes_',phase{1},'_GeneRegulatoryNetwork_Treg_KNOWN_Cell', num2str(CellIndex), '_Treg_FILTERED']);
                plot_biograph(CutCutNewMODfitness,threshold,CUTmRNA_names,FigName,CutCututdcgains);
            end


        end
    end
end
%%%%%%%%%%%%%%%%%%%% Plotting Biograph with different probes piled up in one node %%%%%%%%%%%%%%%%%%%
