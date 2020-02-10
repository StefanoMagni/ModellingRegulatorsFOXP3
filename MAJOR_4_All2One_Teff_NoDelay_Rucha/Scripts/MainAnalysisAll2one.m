%% MAIN SCRIPT TO ANALYSE ALL2(CTLA4,FOXP3,ICOS,IKZF4)
clc;
clear;
close all;

% Select gene of interest to run All2one analysisclear

%--------------
ICOS = 0;  
FOXP3 = 1;
CTLA4 = 0;
IKZF4 = 0;
%--------------
FOXP3Union = 0;
CTLA4Union = 0;
IKZF4Union = 0;
%--------------
FOXP3Window = 0;
CTLA4Window = 0;
IKZF4Window = 0;
%--------------
Treg = 1;
Teff = 0;

Cell_1 = 1;
Cell_2 = 0;
%------------
PlotBioGraph = 1;
PlotDiGraph = 1;

% To exclude genes with mRNA expression level less than the avg. expression level of FOXP3
FOXP3_maxIntesityFilter = 0; 
%--------------



%% Load filtered data of Treg or Teff
if Treg == 1
    CellType = 'Treg';
    load('FILTERED_ExtractedCutData_TregALLgen.mat')
    % Select no. of cut-off genes for the gene of interest
    CutoffGenesNumberFOXP3 = 565;

    CutoffGenesNumberICOS = 540;

    CutoffGenesNumberCTLA4_1 = 450;
    CutoffGenesNumberCTLA4_2 = 500;
    CutoffGenesNumberCTLA4_3 = 400;
    CutoffGenesNumberCTLA4_5 = 550;

    CutoffGenesNumberIKZF4_2 = 259;
    CutoffGenesNumberIKZF4_3 = 470;
    CutoffGenesNumberIKZF4_4 = 600;

    %------------
    CutoffGenesNumberFOXP3_UNION = 500;
    CutoffGenesNumberCTLA4_UNION = 163;
    CutoffGenesNumberIKZF4_UNION = 240;

    %------------
    CutoffGenesNumberFOXP3_2_Window = 23;%20;
    CutoffGenesNumberFOXP3_3_Window = 18;%15;

    CutoffGenesNumberCTLA4_1_Window = 25;
    CutoffGenesNumberCTLA4_2_Window = 20;
    CutoffGenesNumberCTLA4_3_Window = 30;
    CutoffGenesNumberCTLA4_5_Window = 15;

    CutoffGenesNumberIKZF4_2_Window = 25;
    CutoffGenesNumberIKZF4_3_Window = 15;
    CutoffGenesNumberIKZF4_4_Window = 15;
end


%===============%===============%===============
if Teff == 1
    CellType = 'Teff';
    load('FILTERED_ExtractedCutData_TeffALLgen.mat')
    
        % Select no. of cut-off genes for the gene of interest
    CutoffGenesNumberFOXP3 = 1030; %400;

    CutoffGenesNumberICOS = 1370;

    CutoffGenesNumberCTLA4_1 = 1270;
    CutoffGenesNumberCTLA4_2 = 1400;
    CutoffGenesNumberCTLA4_3 = 1200;
    CutoffGenesNumberCTLA4_5 = 1560;

    CutoffGenesNumberIKZF4_2 = 1000;
    CutoffGenesNumberIKZF4_3 = 930;
    CutoffGenesNumberIKZF4_4 = 900;

    %------------
    CutoffGenesNumberFOXP3_UNION = 500;
    CutoffGenesNumberCTLA4_UNION = 163;
    CutoffGenesNumberIKZF4_UNION = 240;

    %------------
    CutoffGenesNumberFOXP3_2_Window = 23;%20;
    CutoffGenesNumberFOXP3_3_Window = 18;%15;

    CutoffGenesNumberCTLA4_1_Window = 25;
    CutoffGenesNumberCTLA4_2_Window = 20;
    CutoffGenesNumberCTLA4_3_Window = 30;
    CutoffGenesNumberCTLA4_5_Window = 15;

    CutoffGenesNumberIKZF4_2_Window = 25;
    CutoffGenesNumberIKZF4_3_Window = 15;
    CutoffGenesNumberIKZF4_4_Window = 15; 
end

%-------------


%------------------------------------------------
% Computing threshold for intensity filter BUT in order to keep FOXP3
load('FOXP3_TregCell12_Probes23.mat')
% maxvalFOXP3_PROBE2 = max(FOXP3_Probe2_Treg12');
avgFOXP3_PROBE2 = mean(FOXP3_Probe2_Treg12');
% maxvalFOXP3_PROBE3 = max(FOXP3_Probe3_Treg12');
avgFOXP3_PROBE3 = mean(FOXP3_Probe3_Treg12');
if FOXP3_maxIntesityFilter == 1
    ThresholdAverageIntensityCell1 = round(min(avgFOXP3_PROBE2(1), avgFOXP3_PROBE3(1)));
    ThresholdAverageIntensityCell2 = round(min(avgFOXP3_PROBE2(2), avgFOXP3_PROBE3(2)));
else
    ThresholdAverageIntensityCell1 = 1;
    ThresholdAverageIntensityCell2 = 1;
end
%------------------------------------------------

%%%%%%%%%% ANALYSIS I %%%%%%%%%%%

if ICOS == 1 
   AnalyseAll2ICOS
end
if FOXP3 == 1
   AnalyseAll2FOXP3
end
if CTLA4 == 1
   AnalyseAll2CTLA4
end
if IKZF4 == 1
   AnalyseAll2IKZF4
end

%%%%%%%%%% ANALYSIS II %%%%%%%%%%%

if FOXP3Union == 1
   AnalyseAll2FOXP3Union
end   
if CTLA4Union == 1
   AnalyseAll2CTLA4Union
end
if IKZF4Union == 1
   AnalyseAll2IKZF4Union
end

%%%%%%%%%% ANALYSIS III %%%%%%%%%%%

if FOXP3Window == 1
   AnalyseAll2FOXP3Window
end   
if CTLA4Window == 1
   AnalyseAll2CTLA4Window
end
if IKZF4Window == 1
   AnalyseAll2IKZF4Window
end

