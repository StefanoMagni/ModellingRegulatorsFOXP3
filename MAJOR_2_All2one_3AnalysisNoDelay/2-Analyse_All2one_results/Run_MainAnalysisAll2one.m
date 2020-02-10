%% MAIN SCRIPT TO ANALYSE ALL2(CTLA4,FOXP3,ICOS,IKZF4)
clc;
clear;
close all;

% Select gene of interest to run All2one analysis

%--------------
ICOS = 1;
FOXP3 = 1;
CTLA4 = 1;
IKZF4 = 1;
%--------------
FOXP3Union = 1;
CTLA4Union = 1;
IKZF4Union = 1;
%--------------
FOXP3Window = 1;
CTLA4Window = 1;
IKZF4Window = 1;
%--------------

Treg = 1;
Teff = 0;

if Treg == 1
    CellType = 'Treg';
end

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
%------------
PlotBioGraph = 1;

%------------------------------------------------
% Computing threshold for intensity filter BUT in order to keep FOXP3
load('FOXP3_TregCell12_Probes23.mat')
% maxvalFOXP3_PROBE2 = max(FOXP3_Probe2_Treg12');
avgFOXP3_PROBE2 = mean(FOXP3_Probe2_Treg12');
% maxvalFOXP3_PROBE3 = max(FOXP3_Probe3_Treg12');
avgFOXP3_PROBE3 = mean(FOXP3_Probe3_Treg12');
ThresholdAverageIntensityCell1 = round(min(avgFOXP3_PROBE2(1), avgFOXP3_PROBE3(1)));
ThresholdAverageIntensityCell2 = round(min(avgFOXP3_PROBE2(2), avgFOXP3_PROBE3(2)));
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

