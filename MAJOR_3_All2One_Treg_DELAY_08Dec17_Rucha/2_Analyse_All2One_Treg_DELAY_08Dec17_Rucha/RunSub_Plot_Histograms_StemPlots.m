%% PLOT HISTOGRAMS
% To display distribution of fitness score against the no. of genes
% for each of the 6 delay cases
% 11 Dec, 2017 (By Rucha Sawlekar)


    load('FitnessFor_0Delay.mat')
    MaxVal = max(FitnessFor_0Delay);
    figure
    subplot(1,2,1)
    stem(FitnessFor_0Delay,'DisplayName','FitnessFor_0Delay')
    xlim([0 3515]); xlabel('Genes'); ylabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',11); 
         title('Fitness for 0 Delay', 'FontWeight','bold', 'FontSize',13)
    subplot(1,2,2)
    histogram(FitnessFor_0Delay(:,1),20,'FaceAlpha',0.5); hold on;
    histogram(FitnessFor_0Delay(:,2),20); hold on;
    histogram(FitnessFor_0Delay(:,3),20); hold on;
    histogram(FitnessFor_0Delay(:,4),20); hold on;
    ylabel('No. of Genes'); xlabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',11); 
%     labels = arrayfun(@(value) num2str(value,'%2.2f'),cell2mat(FitnessFor_0Delay(:,1:4)),'UniformOutput',false);
%     text(Delay,cell2mat(FitnessFor_0Delay(:,1:4)),labels,'HorizontalAlignment','center','VerticalAlignment','bottom')    

    title('Fitness for 0 Delay', 'FontWeight','bold', 'FontSize',15)

    %-------------

    load('FitnessFor_20Delay.mat')
    MaxVal = max(FitnessFor_20Delay);
    figure
    subplot(1,2,1)
    stem(FitnessFor_20Delay,'DisplayName','FitnessFor_20Delay')
    xlim([0 3515]); xlabel('Genes'); ylabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',14); 
         title('Fitness for 20 Delay', 'FontWeight','bold', 'FontSize',13)
    subplot(1,2,2)
    histogram(FitnessFor_20Delay(:,1),20,'FaceAlpha',0.5); hold on;
    histogram(FitnessFor_20Delay(:,2),20); hold on;
    histogram(FitnessFor_20Delay(:,3),20); hold on;
    histogram(FitnessFor_20Delay(:,4),20); hold on;
    ylabel('No. of Genes'); xlabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',14); 
    title('Fitness for 20 Delay', 'FontWeight','bold', 'FontSize',15)

    %-------------

    load('FitnessFor_40Delay.mat')
    MaxVal = max(FitnessFor_40Delay);
    figure
    subplot(1,2,1)
    stem(FitnessFor_40Delay,'DisplayName','FitnessFor_40Delay')
    xlim([0 3515]); xlabel('Genes'); ylabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',14); 
         title('Fitness for 40 Delay', 'FontWeight','bold', 'FontSize',15)
    subplot(1,2,2)
    histogram(FitnessFor_40Delay(:,1),20,'FaceAlpha',0.5); hold on;
    histogram(FitnessFor_40Delay(:,2),20); hold on;
    histogram(FitnessFor_40Delay(:,3),20); hold on;
    histogram(FitnessFor_40Delay(:,4),20); hold on;
    ylabel('No. of Genes'); xlabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',14); 
    title('Fitness for 40 Delay', 'FontWeight','bold', 'FontSize',15)

    %-------------

    load('FitnessFor_60Delay.mat')
    MaxVal = max(FitnessFor_60Delay);
    figure
    subplot(1,2,1)
    stem(FitnessFor_60Delay,'DisplayName','FitnessFor_60Delay')
    xlim([0 3515]); xlabel('Genes'); ylabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',14); 
         title('Fitness for 60 Delay', 'FontWeight','bold', 'FontSize',15)
    subplot(1,2,2)
    histogram(FitnessFor_60Delay(:,1),20,'FaceAlpha',0.5); hold on;
    histogram(FitnessFor_60Delay(:,2),20); hold on;
    histogram(FitnessFor_60Delay(:,3),20); hold on;
    histogram(FitnessFor_60Delay(:,4),20); hold on;
    ylabel('No. of Genes'); xlabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',14); 
    title('Fitness for 60 Delay', 'FontWeight','bold', 'FontSize',15)

    %-------------

    load('FitnessFor_80Delay.mat')
    MaxVal = max(FitnessFor_80Delay);
    figure
    subplot(1,2,1)
    stem(FitnessFor_80Delay,'DisplayName','FitnessFor_80Delay')
    xlim([0 3515]); xlabel('Genes'); ylabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',14); 
         title('Fitness for 80 Delay', 'FontWeight','bold', 'FontSize',15)
    subplot(1,2,2)
    histogram(FitnessFor_80Delay(:,1),20,'FaceAlpha',0.5); hold on;
    histogram(FitnessFor_80Delay(:,2),20); hold on;
    histogram(FitnessFor_80Delay(:,3),20); hold on;
    histogram(FitnessFor_80Delay(:,4),20); hold on;
    ylabel('No. of Genes'); xlabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',14); 
    title('Fitness for 80 Delay', 'FontWeight','bold', 'FontSize',15)

    %------------

    load('FitnessFor_100Delay.mat')
    MaxVal = max(FitnessFor_100Delay);
    figure
    subplot(1,2,1)
    stem(FitnessFor_100Delay,'DisplayName','FitnessFor_100Delay')
    xlim([0 3515]); xlabel('Genes'); ylabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',14); 
         title('Fitness for 100 Delay', 'FontWeight','bold', 'FontSize',15)
    subplot(1,2,2)
    histogram(FitnessFor_100Delay(:,1),20,'FaceAlpha',0.5); hold on;
    histogram(FitnessFor_100Delay(:,2),20); hold on;
    histogram(FitnessFor_100Delay(:,3),20); hold on;
    histogram(FitnessFor_100Delay(:,4),20); hold on;
    ylabel('No. of Genes'); xlabel('Fitness Score')
    legend({['P2-C1, MaxFit=' num2str(MaxVal(1))],['P3-C1, MaxFit=' num2str(MaxVal(2))],...
            ['P2-C2, MaxFit=' num2str(MaxVal(3))],['P3-C2, MaxFit=' num2str(MaxVal(4))]},...
             'FontWeight','bold', 'FontSize',14); 
    title('Fitness for 100 Delay', 'FontWeight','bold', 'FontSize',15)
    
    
    
    
%       ******************           END            ******************
    


    