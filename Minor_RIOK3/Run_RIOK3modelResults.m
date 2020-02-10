% Look up RIOK3 (Probe ID = 202130_at) model results
% It pops up as a potential upstream regulator (repressor) of FOXP3 probe 2
% and 3

% CAREFULL !!!!!!!
% THIS SCRIPT MAY CREATE 1 FINAL FIGURE OF SZE ABOVE 1 GB FOR UNEXPLAINED
% RESONS...

clear; close all; clc;
warning off;

%%%%%%%%%%% FOXP3 Probe 2 %%%%%%%%%%

%%%%% Cell 1

load('FOXP3_Results_Treg_AllGenes_Cell1ProbeSet_2.mat', 'models');
modelsProbe2Cell1 = models;
load('FOXP3_Results_Treg_AllGenes_Cell1ProbeSet_2.mat', 'IOs');
IOsProbe2Cell1 = IOs;
load('FILTERED_ExtractedCutData_TregALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell1');
TableOfExtractedGeneNameAndIDs_Probe2Cell1 = TableOfExtractedGeneNameAndIDs_cell1;

temp = strcmp('202130_at',TableOfExtractedGeneNameAndIDs_Probe2Cell1(:,1));
idx = find(temp == 1)
modelsProbe2Cell1{idx}

hl=figure; bode(modelsProbe2Cell1{idx}); saveas(hl,[pwd '/BodeP2C1.fig']);
hl=figure; step(modelsProbe2Cell1{idx}); saveas(hl,[pwd '/StepP2C1.fig']);
hl=figure; compare(IOsProbe2Cell1{idx},modelsProbe2Cell1{idx}); saveas(hl,[pwd '/CompareP2C1.fig']);
hl=figure; plot(IOsProbe2Cell1{idx}); saveas(hl,[pwd '/IOsP2C1.fig']);

%%%%% Cell 2

load('FOXP3_Results_Treg_AllGenes_Cell2ProbeSet_2.mat', 'models');
modelsProbe2Cell2 = models;
load('FOXP3_Results_Treg_AllGenes_Cell2ProbeSet_2.mat', 'IOs');
IOsProbe2Cell2 = IOs;
load('FILTERED_ExtractedCutData_TregALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell2');
TableOfExtractedGeneNameAndIDs_Probe2Cell2 = TableOfExtractedGeneNameAndIDs_cell2;

temp = strcmp('202130_at',TableOfExtractedGeneNameAndIDs_Probe2Cell2(:,1));
idx = find(temp == 1)
modelsProbe2Cell2{idx}

hl=figure;bode(modelsProbe2Cell2{idx}); saveas(hl,[pwd '/BodeP2C2.fig']);
hl=figure;step(modelsProbe2Cell2{idx}); saveas(hl,[pwd '/StepP2C2.fig']);
hl=figure;compare(IOsProbe2Cell2{idx},modelsProbe2Cell2{idx}); saveas(hl,[pwd '/CompareP2C2.fig']);
hl=figure;plot(IOsProbe2Cell2{idx}); saveas(hl,[pwd '/IOsP2C2.fig']);

%%%%%%%%%%% FOXP3 Probe 2 %%%%%%%%%%

%%%%% Cell 1

load('FOXP3_Results_Treg_AllGenes_Cell1ProbeSet_3.mat', 'models');
modelsProbe3Cell1 = models;
load('FOXP3_Results_Treg_AllGenes_Cell1ProbeSet_3.mat', 'IOs');
IOsProbe3Cell1 = IOs;
load('FILTERED_ExtractedCutData_TregALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell1');
TableOfExtractedGeneNameAndIDs_Probe3Cell1 = TableOfExtractedGeneNameAndIDs_cell1;

temp = strcmp('202130_at',TableOfExtractedGeneNameAndIDs_Probe3Cell1(:,1));
idx = find(temp == 1)
modelsProbe3Cell1{idx}

hl=figure;bode(modelsProbe3Cell1{idx}); saveas(hl,[pwd '/BodeP3C1.fig']);
hl=figure;step(modelsProbe3Cell1{idx}); saveas(hl,[pwd '/StepP3C1.fig']);
hl=figure;compare(IOsProbe3Cell1{idx},modelsProbe3Cell1{idx}); saveas(hl,[pwd 'CompareP3C1.fig']);
hl=figure;plot(IOsProbe3Cell1{idx}); saveas(hl,[pwd '/IOsP3C1.fig']);

%%%%% Cell 2

load('FOXP3_Results_Treg_AllGenes_Cell2ProbeSet_3.mat', 'models');
modelsProbe3Cell2 = models;
load('FOXP3_Results_Treg_AllGenes_Cell2ProbeSet_3.mat', 'IOs');
IOsProbe3Cell2 = IOs;
load('FILTERED_ExtractedCutData_TregALLgen.mat', 'TableOfExtractedGeneNameAndIDs_cell2');
TableOfExtractedGeneNameAndIDs_Probe3Cell2 = TableOfExtractedGeneNameAndIDs_cell2;

temp = strcmp('202130_at',TableOfExtractedGeneNameAndIDs_Probe3Cell2(:,1));
idx = find(temp == 1)
modelsProbe3Cell2{idx}

hl=figure;bode(modelsProbe3Cell2{idx}); saveas(hl,[pwd '/BodeP3C2.fig']);
hl=figure;step(modelsProbe3Cell2{idx}); saveas(hl,[pwd '/StepP3C2.fig']);
hl=figure;compare(IOsProbe3Cell2{idx},modelsProbe3Cell2{idx}); saveas(hl,[pwd '/CompareP3C2.fig']);
hl=figure;plot(IOsProbe3Cell2{idx}); saveas(hl,[pwd '/IOsP3C2.fig']);

warning on;



%%%%%%%%%%%%%%%%%%%%%%% Extract Other Models for comparison %%%%%%%%%%%%%%%%%%%%%%%%%
%%% FOXP3 Probe 2 - LOOK FOR IKZF4 as INPUT %%%
%%%%% Cell 1
temp = strcmp('226759_at',TableOfExtractedGeneNameAndIDs_Probe2Cell1(:,1));
idx = find(temp == 1)
modelsProbe2Cell1{idx}
%%%%% Cell 2
temp = strcmp('229752_at',TableOfExtractedGeneNameAndIDs_Probe2Cell2(:,1));
idx = find(temp == 1)
modelsProbe2Cell2{idx}
%%%%% Cell 2
temp = strcmp('226761_at',TableOfExtractedGeneNameAndIDs_Probe2Cell2(:,1));
idx = find(temp == 1)
modelsProbe2Cell2{idx}
% Define transfer functions
%%% END %%%
%%% FOXP3 probes 2/3 look at SARNP as input %%%
%%%%% Cell 1
temp = strcmp('224914_s_at',TableOfExtractedGeneNameAndIDs_Probe2Cell1(:,1));
idx = find(temp == 1)
modelsProbe2Cell1{idx}
%%%%% Cell 2
temp = strcmp('224914_s_at',TableOfExtractedGeneNameAndIDs_Probe2Cell2(:,1));
idx = find(temp == 1)
modelsProbe2Cell2{idx}
%%%%% Cell 1
temp = strcmp('224914_s_at',TableOfExtractedGeneNameAndIDs_Probe3Cell1(:,1));
idx = find(temp == 1)
modelsProbe3Cell1{idx}
%%%%% Cell 2
temp = strcmp('224914_s_at',TableOfExtractedGeneNameAndIDs_Probe3Cell2(:,1));
idx = find(temp == 1)
modelsProbe3Cell2{idx}
%%% END %%%



%%%%%%%%%%%%%%%%%%%%%%% Compute vgap %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Probe_2 - Cell_1 VS Cell_2 %%%%%%%%%%

%How to use vgap function:
%H1p2c1(s) = -1.25/(s + 0.9631); H2p2c2(s) = -1.273/(s + 1.616);
% Define transfer functions
% RIOK3
N1 = {[-1.25]}; D1 = {[1  0.9631]}; Hprobe2cell1 = tf(N1,D1);
N2 = {[-1.273]}; D2 = {[1 1.616]}; Hprobe2cell2 = tf(N2,D2); 
N3 = {[-1.509]}; D3 = {[1  1.411]}; Hprobe3cell1 = tf(N3,D3);
N4 = {[-1.138]}; D4 = {[1 1.776]}; Hprobe3cell2 = tf(N4,D4); 

% IKZF4
N1 = {[0.6712]}; D1 = {[1  0.09093]}; Hprobe2cell1inputIKZF4 = tf(N1,D1);
N2 = {[0.9074]}; D2 = {[1 0.7088]}; Hprobe2cell2inputIKZF4a = tf(N2,D2); 
N3 = {[1.532]}; D3 = {[1  1.949]}; Hprobe2cell2inputIKZF4b = tf(N3,D3);

%Compute gaps + frequencies
%Range of frequencies (give it large, so that you can see the whole
%behaviour)
f1 = 10 % 100; %2; %Upper bound (Hz)
f2 = 1/10 %1/100; %1/20; %Lower bound (Hz)

% Compute nu gap
[nugapX1, freqX1, distnuX1, wX1]=richgap(Hprobe2cell1inputIKZF4, Hprobe2cell2inputIKZF4a, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapX2, freqX2, distnuX2, wX2]=richgap(Hprobe2cell1inputIKZF4, Hprobe2cell2inputIKZF4b, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapX3, freqX3, distnuX3, wX3]=richgap(Hprobe2cell2inputIKZF4a, Hprobe2cell2inputIKZF4b, 2*pi/(1/f1), 2*pi/(1/f2));

[nugapZ1, freqZ1, distnuZ1, wZ1]=richgap(Hprobe2cell1inputIKZF4, Hprobe2cell1, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapZ2, freqZ2, distnuZ2, wZ2]=richgap(Hprobe2cell2inputIKZF4a, Hprobe2cell2, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapZ3, freqZ3, distnuZ3, wZ3]=richgap(Hprobe2cell2inputIKZF4b, Hprobe2cell2, 2*pi/(1/f1), 2*pi/(1/f2));

[nugapZ4, freqZ4, distnuZ4, wZ4]=richgap(Hprobe2cell1inputIKZF4, Hprobe2cell2, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapZ5, freqZ5, distnuZ5, wZ5]=richgap(Hprobe2cell2inputIKZF4a, Hprobe2cell1, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapZ6, freqZ6, distnuZ6, wZ6]=richgap(Hprobe2cell2inputIKZF4b, Hprobe2cell1, 2*pi/(1/f1), 2*pi/(1/f2));

[nugapZ7, freqZ7, distnuZ7, wZ7]=richgap(Hprobe2cell1inputIKZF4, Hprobe3cell1, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapZ8, freqZ8, distnuZ8, wZ8]=richgap(Hprobe2cell2inputIKZF4a, Hprobe3cell2, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapZ9, freqZ9, distnuZ9, wZ9]=richgap(Hprobe2cell2inputIKZF4b, Hprobe3cell2, 2*pi/(1/f1), 2*pi/(1/f2));

[nugapZ10, freqZ10, distnuZ10, wZ10]=richgap(Hprobe2cell1inputIKZF4, Hprobe3cell2, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapZ11, freqZ11, distnuZ11, wZ11]=richgap(Hprobe2cell2inputIKZF4a, Hprobe3cell1, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapZ12, freqZ12, distnuZ12, wZ12]=richgap(Hprobe2cell2inputIKZF4b, Hprobe3cell1, 2*pi/(1/f1), 2*pi/(1/f2));

% SARNP
N1 = {[-2.037]}; D1 = {[1  2.575]}; Hprobe2cell1inputSARNP = tf(N1,D1);
N2 = {[-0.9412]}; D2 = {[1 0.9697]}; Hprobe2cell2inputSARNP = tf(N2,D2); 
N3 = {[-5.782]}; D3 = {[1  8.382]}; Hprobe3cell1inputSARNP = tf(N3,D3);
N4 = {[-3.129]}; D4 = {[1 4.294]}; Hprobe3cell2inputSARNP = tf(N4,D4);



hl = figure;

[nugap, freq, distnu, w]=richgap(Hprobe2cell1, Hprobe2cell2, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapSARNP, freqSARNP, distnuSARNP, wSARNP]=richgap(Hprobe2cell1inputSARNP, Hprobe2cell2inputSARNP, 2*pi/(1/f1), 2*pi/(1/f2));

fprintf('nugap = %d at freq = %d \n',nugap, freq/(2*pi));
subplot(3,2,1); 
a = plot(w./(2*pi),distnu,'LineWidth',3);  hold on;
b = plot(wX1./(2*pi),distnuX1,':','LineWidth',1.5,'color', 'green'); hold on; ...
    c = plot(wX2./(2*pi),distnuX2,'--','color', 'green'); hold on; ...
    d = plot(wX3./(2*pi),distnuX3,'color', 'green'); hold on;
e = plot(w./(2*pi),distnuZ1,'color', 'red');  hold on; plot(w./(2*pi),distnuZ2,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ3,'color', 'red');  hold on; plot(w./(2*pi),distnuZ4,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ5,'color', 'red');  hold on; plot(w./(2*pi),distnuZ6,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ7,'color', 'red');  hold on; plot(w./(2*pi),distnuZ8,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ9,'color', 'red');  hold on; plot(w./(2*pi),distnuZ10,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ11,'color', 'red');  hold on; plot(w./(2*pi),distnuZ12,'color', 'red');
f = plot(w./(2*pi),distnuSARNP,'color', 'black');  hold on;
xlabel('Freq. (Hz)'); ylabel('nugaps'); grid; grid minor; box on;
set(gca, 'XScale', 'log');
xlim([0.01 100]);
ylim([0.0 0.5]);
title('Cell 1 VS Cell 2 (probe 2)');
legend([a b c d e f], 'RIOK3             ', 'IKZF4 PaC1 VS PbC2', 'IKZF4 PaC1 VS PcC2', ...
    'IKZF4 PbC2 VS PcC2', 'IKZF4 VS RIOK3 var', 'SARNP             ')

[nugap, freq, distnu, w]=richgap(Hprobe3cell1, Hprobe3cell2, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapSARNP, freqSARNP, distnuSARNP, wSARNP]=richgap(Hprobe3cell1inputSARNP, Hprobe3cell2inputSARNP, 2*pi/(1/f1), 2*pi/(1/f2));

fprintf('nugap = %d at freq = %d \n',nugap, freq/(2*pi));
subplot(3,2,2);
plot(w./(2*pi),distnu,'LineWidth',3);  hold on;
plot(wX1./(2*pi),distnuX1,':','LineWidth',1.5,'color', 'green'); hold on; plot(wX2./(2*pi),distnuX2,'--','color', 'green'); ...
    hold on; plot(wX3./(2*pi),distnuX3,'color', 'green');
plot(w./(2*pi),distnuZ1,'color', 'red');  hold on; plot(w./(2*pi),distnuZ2,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ3,'color', 'red');  hold on; plot(w./(2*pi),distnuZ4,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ5,'color', 'red');  hold on; plot(w./(2*pi),distnuZ6,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ7,'color', 'red');  hold on; plot(w./(2*pi),distnuZ8,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ9,'color', 'red');  hold on; plot(w./(2*pi),distnuZ10,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ11,'color', 'red');  hold on; plot(w./(2*pi),distnuZ12,'color', 'red');
plot(w./(2*pi),distnuSARNP,'color', 'black');  hold on;
xlabel('Freq. (Hz)'); ylabel('nugaps'); grid; grid minor; box on;
set(gca, 'XScale', 'log');
xlim([0.01 100]);
ylim([0.0 0.5]);
title('Cell 1 VS Cell 2 (probe 3)');

[nugap, freq, distnu, w]=richgap(Hprobe2cell1, Hprobe3cell1, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapSARNP, freqSARNP, distnuSARNP, wSARNP]=richgap(Hprobe2cell1inputSARNP, Hprobe3cell2inputSARNP, 2*pi/(1/f1), 2*pi/(1/f2));

fprintf('nugap = %d at freq = %d \n',nugap, freq/(2*pi));
subplot(3,2,3);
plot(w./(2*pi),distnu,'LineWidth',3);  hold on;
plot(wX1./(2*pi),distnuX1,':','LineWidth',1.5,'color', 'green'); hold on; plot(wX2./(2*pi),distnuX2,'--','color', 'green'); ...
    hold on; plot(wX3./(2*pi),distnuX3,'color', 'green');
plot(w./(2*pi),distnuZ1,'color', 'red');  hold on; plot(w./(2*pi),distnuZ2,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ3,'color', 'red');  hold on; plot(w./(2*pi),distnuZ4,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ5,'color', 'red');  hold on; plot(w./(2*pi),distnuZ6,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ7,'color', 'red');  hold on; plot(w./(2*pi),distnuZ8,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ9,'color', 'red');  hold on; plot(w./(2*pi),distnuZ10,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ11,'color', 'red');  hold on; plot(w./(2*pi),distnuZ12,'color', 'red');
plot(w./(2*pi),distnuSARNP,'color', 'black');  hold on;
xlabel('Freq. (Hz)'); ylabel('nugaps'); grid; grid minor; box on;
set(gca, 'XScale', 'log');
xlim([0.01 100]);
ylim([0.0 0.5]);
title('Probe 2 VS Probe 3 (Cell 1)');

[nugap, freq, distnu, w]=richgap(Hprobe2cell2, Hprobe3cell2, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapSARNP, freqSARNP, distnuSARNP, wSARNP]=richgap(Hprobe2cell2inputSARNP, Hprobe3cell2inputSARNP, 2*pi/(1/f1), 2*pi/(1/f2));

fprintf('nugap = %d at freq = %d \n',nugap, freq/(2*pi));
subplot(3,2,4);
plot(w./(2*pi),distnu,'LineWidth',3);  hold on;
plot(wX1./(2*pi),distnuX1,':','LineWidth',1.5,'color', 'green'); hold on; plot(wX2./(2*pi),distnuX2,'--','color', 'green'); ...
    hold on; plot(wX3./(2*pi),distnuX3,'color', 'green');
plot(w./(2*pi),distnuZ1,'color', 'red');  hold on; plot(w./(2*pi),distnuZ2,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ3,'color', 'red');  hold on; plot(w./(2*pi),distnuZ4,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ5,'color', 'red');  hold on; plot(w./(2*pi),distnuZ6,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ7,'color', 'red');  hold on; plot(w./(2*pi),distnuZ8,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ9,'color', 'red');  hold on; plot(w./(2*pi),distnuZ10,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ11,'color', 'red');  hold on; plot(w./(2*pi),distnuZ12,'color', 'red');
plot(w./(2*pi),distnuSARNP,'color', 'black');  hold on;
xlabel('Freq. (Hz)'); ylabel('nugaps'); grid; grid minor; box on;
set(gca, 'XScale', 'log');
xlim([0.01 100]);
ylim([0.0 0.5]);
title('Probe 2 VS Probe 3 (Cell 2)');

[nugap, freq, distnu, w]=richgap(Hprobe2cell1, Hprobe3cell2, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapSARNP, freqSARNP, distnuSARNP, wSARNP]=richgap(Hprobe2cell1inputSARNP, Hprobe3cell2inputSARNP, 2*pi/(1/f1), 2*pi/(1/f2));

fprintf('nugap = %d at freq = %d \n',nugap, freq/(2*pi));
subplot(3,2,5);
plot(w./(2*pi),distnu,'LineWidth',3);  hold on;
plot(wX1./(2*pi),distnuX1,':','LineWidth',1.5,'color', 'green'); hold on; plot(wX2./(2*pi),distnuX2,'--','color', 'green'); ...
    hold on; plot(wX3./(2*pi),distnuX3,'color', 'green');
plot(w./(2*pi),distnuZ1,'color', 'red');  hold on; plot(w./(2*pi),distnuZ2,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ3,'color', 'red');  hold on; plot(w./(2*pi),distnuZ4,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ5,'color', 'red');  hold on; plot(w./(2*pi),distnuZ6,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ7,'color', 'red');  hold on; plot(w./(2*pi),distnuZ8,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ9,'color', 'red');  hold on; plot(w./(2*pi),distnuZ10,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ11,'color', 'red');  hold on; plot(w./(2*pi),distnuZ12,'color', 'red');
plot(w./(2*pi),distnuSARNP,'color', 'black');  hold on;
xlabel('Freq. (Hz)'); ylabel('nugaps'); grid; grid minor; box on;
set(gca, 'XScale', 'log');
xlim([0.01 100]);
ylim([0.0 0.5]);
title('Probe 2 Cell 1 VS Probe 3 Cell 2');

[nugap, freq, distnu, w]=richgap(Hprobe2cell2, Hprobe3cell1, 2*pi/(1/f1), 2*pi/(1/f2));
[nugapSARNP, freqSARNP, distnuSARNP, wSARNP]=richgap(Hprobe2cell2inputSARNP, Hprobe3cell1inputSARNP, 2*pi/(1/f1), 2*pi/(1/f2));

fprintf('nugap = %d at freq = %d \n',nugap, freq/(2*pi));
subplot(3,2,6);
plot(w./(2*pi),distnu,'LineWidth',3);  hold on;
plot(wX1./(2*pi),distnuX1,':','LineWidth',1.5,'color', 'green'); hold on; plot(wX2./(2*pi),distnuX2,'--','color', 'green'); ...
    hold on; plot(wX3./(2*pi),distnuX3,'color', 'green');
plot(w./(2*pi),distnuZ1,'color', 'red');  hold on; plot(w./(2*pi),distnuZ2,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ3,'color', 'red');  hold on; plot(w./(2*pi),distnuZ4,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ5,'color', 'red');  hold on; plot(w./(2*pi),distnuZ6,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ7,'color', 'red');  hold on; plot(w./(2*pi),distnuZ8,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ9,'color', 'red');  hold on; plot(w./(2*pi),distnuZ10,'color', 'red');  hold on;
plot(w./(2*pi),distnuZ11,'color', 'red');  hold on; plot(w./(2*pi),distnuZ12,'color', 'red');
plot(w./(2*pi),distnuSARNP,'color', 'black');  hold on;
xlabel('Freq. (Hz)'); ylabel('nugaps'); grid; grid minor; box on;
set(gca, 'XScale', 'log');
xlim([0.01 100]);
ylim([0.0 0.5]);
title('Probe 2 Cell 2 VS Probe 3 Cell 1');

suptitle('Gaps VS Frequencies');

saveas(hl,[pwd '/RIOK3results/NuGapAnalysis.fig']);


