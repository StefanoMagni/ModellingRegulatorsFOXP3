
#### FOLDER == ‘MAJOR_3_All2One_Treg_DELAY_08Dec17_Rucha’
The following information is regarding folder mentioned above in main folder on owncloud, ‘TregProjectSCG_MASTER_FOLDER’. Data to run these codes and results produced by these codes are stored in the same file name, ‘MAJOR_3_All2One_Treg_DELAY_08Dec17_Rucha’, on owncloud.
Following codes written by RUCHA SAWLEKAR (Dec2017-Jan2018)Sequence to follow:1.	To first intersect then take union of Cell-1 and 2. This results in 3515 genes each for P2C1, P3C1, P2C2, P3C2Code to RUN => UnionOfCell1and2_IMP_27Nov.m2.	To summarize and save the information on fitness score after all2one of 6 delay casesCode to RUN => Seach_All_Cases.m 3.	Search and plot info in all 6-delay cases (Treg-1&2) for manual selectionCode to RUN => Seach_All_Cases_PLOT.m4.	Now we grade the genes as ‘1’, ‘1.5’, ‘2’ etc. and now plot only those for selected grade, for eg: plot only for ‘Grade=1’ or only for ‘Grade=2’Code to RUN => Plot_GradedGenes.m (Also IMP code: RankFitness_DELAY_HISTG_11Dec.m)
Main mat files to refer for given grades are:
- GradedGenes_20Jan18.mat   and
- GradingGenes_07Dec.mat5.	To plot only ‘---’ gene namesCode to RUN => PlotOnlyDashGenes.m (Also IMP code: RankFitness_DELAY_HISTG_11Dec.m)
