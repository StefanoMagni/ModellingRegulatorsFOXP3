%% MAIN SCRIPT TO ANALYSE ALL2(CTLA4,FOXP3,ICOS,IKZF4)
clc;
clear;
close all;

% Select gene of interest to run All2one analysisclear

%--------------
FOXP3 = 1;

%--------------
Treg = 1;
Teff = 0;

Cell_1 = 0;
Cell_2 = 1;
%--------------

PlotBioGraph = 0;
PlotDiGraph = 0;

SelectCandidates = 2;

%--------------

% Load filtered data of Treg or Teff
if Treg == 1
    CellType = 'Treg';
    FileNameToLoad = ['FILTERED_ExtractedCutData_' CellType 'ALLgen.mat'];        
    load(FileNameToLoad)    
    % Select no. of cut-off genes for the gene of interest
    CutoffGenesNumberFOXP3C1 = size(FILTERED_AFFI_FLAGS_Tcell1,1); %1500;%565;
    CutoffGenesNumberFOXP3C2 = size(FILTERED_AFFI_FLAGS_Tcell2,1);
    
% %     CutoffGenesNumberFOXP3C1 = 1000;
% %     CutoffGenesNumberFOXP3C2 = 1000;

% %     CutoffGenesNumberFOXP3 = 565;    
end

%--------------%--------------%--------------
if Teff == 1
    CellType = 'Teff';
    FileNameToLoad = ['FILTERED_ExtractedCutData_' CellType 'ALLgen.mat'];        
    load(FileNameToLoad)    
    % Select no. of cut-off genes for the gene of interest
    CutoffGenesNumberFOXP3 = 1030; %400;
end


%------------------------------------------------
% Computing threshold for intensity filter BUT in order to keep FOXP3
load('FOXP3_TregCell12_Probes23.mat')
% maxvalFOXP3_PROBE2 = max(FOXP3_Probe2_Treg12');
avgFOXP3_PROBE2 = mean(FOXP3_Probe2_Treg12');
% maxvalFOXP3_PROBE3 = max(FOXP3_Probe3_Treg12');
avgFOXP3_PROBE3 = mean(FOXP3_Probe3_Treg12');
ThresholdAverageIntensityCell1 = 1; % round(min(avgFOXP3_PROBE2(1), avgFOXP3_PROBE3(1)));
ThresholdAverageIntensityCell2 = 1; % round(min(avgFOXP3_PROBE2(2), avgFOXP3_PROBE3(2)));
% By setting threhold = 1, we practically remove " intensity < than average
% of FOXP3" in both cells
%------------------------------------------------

%%%%%%%%%% ANALYSIS I - RANKING %%%%%%%%%%%

if FOXP3 == 1
%    AnalyseAll2FOXP3
   AnalyseAll2FOXP3_Ranking_V2
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      END      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






