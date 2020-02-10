% Main script to scroll across data of Treg/Teff
% Created by Stefano Magni in April 2017, stefano.magni@uni.lu

% THIS SCRIPT DOES NOT SEEMS TO WORK WITH CURRENT VCERSION OF .mat FILES OF
% EXTRACTED GENES DUE TO DIFFERENT AFFI FLAGS.

fromRIGHT = 0;
CLthrasholdRIGHTtail = (0.50);%1-0.005/54000.);

fromLEFT = 1;
CLthrasholdLEFTtail = 0.99;

PassaBanda = 0;
MinCLband = 0.45;
MaxCLband = 0.55;

STEP = 25;

% Feng Filter
FilterMaxAmplitude = 100;
FilterMean = 50;

% Affimetric flags filter
FilterAFFIflagsAbsent = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('ListAFFIMETRICflags.mat')
load('ExtractedCutData_TeffALLgen.mat');
%load('ExtractedCutData_TregALLgen.mat')

QCL = LjungBoxTestQandCL_Tcell1;%[ 10, 0.005; 120, 0.995; 60, 0.54] %
% LjungBoxTestQandCL_Tcell2
TimePoints = Tcell1DataTable;%[ 1, 1, 5, 1, 1; 2, 2, 2, 4, 4; 0, 3, 2, 4, 1];%
% Tcell2DataTable
PositionForAFFIflag = 1; %   2     3     4
                  % Teff1 Teff2 Treg1 Treg2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Names = TableOfExtractedGeneNameAndIDs(:,4);%{'Gene1';'Gene2';'Gene3'}%
ListAFFIMETRICflags = NewListOfAllFLAGSrightOrder;

QCLplusTimePoints = [QCL,TimePoints];
[~, Index] = sort(QCLplusTimePoints(:,1));
SORTED_QCLplusTimePoints = QCLplusTimePoints(Index,:);
SORTED_Names = Names(Index,:);
SORTED_ListAFFIMETRICflags = ListAFFIMETRICflags(Index,:);

Times = [0:20:360];

if fromLEFT == 1
    STOPindex = 0;
    for i = 1 : numel(SORTED_QCLplusTimePoints(:,1))
        if SORTED_QCLplusTimePoints(i,2) <= CLthrasholdLEFTtail
            STOPindex = i;
        else 
            break;
        end
    end
    disp(['Number of selected genes in tail is ' num2str(STOPindex)])
    disp(' ')
    for k= 1 : STEP : STOPindex
        %k = STOPindex - j + 1;
        Q = SORTED_QCLplusTimePoints(k,1);
        CL = SORTED_QCLplusTimePoints(k,2);
        SORTED_Names(k);
        MyTitle = [strjoin(SORTED_Names(k)) ...
            '        Q = ' num2str(Q) ...
            '        CL = ' num2str(CL)];

        disp(['Gene ' num2str(k) ' out of ' num2str(STOPindex)])
        if max(SORTED_QCLplusTimePoints(k,1+2:19+2)) >= FilterMaxAmplitude & ...
                mean(SORTED_QCLplusTimePoints(k,1+2:19+2)) >= FilterMean
            AFFIMETRICflag = SORTED_ListAFFIMETRICflags(k,PositionForAFFIflag);
            disp(AFFIMETRICflag)
            if not(FilterAFFIflagsAbsent) | ...
                    not(strcmp(AFFIMETRICflag{1},'YesABSENT'))
                MyFigure = figure;
                subplot(2,1,1)
                plot(Times,SORTED_QCLplusTimePoints(k,1+2:19+2));
                xlabel('time (min)'); ylabel('Gene Expression');
                title(MyTitle);
                subplot(2,1,2)
                autocorr(SORTED_QCLplusTimePoints(k,1+2:19+2))
                waitfor(MyFigure);
                disp(MyTitle)
                disp(' ')
            end
        else
            disp(['Excluded because of Feng filter'])
            disp(' ')
        end
    end
elseif fromRIGHT == 1
    STOPindex = 0;
    MAXindex = numel(SORTED_QCLplusTimePoints(:,1));
    for i = 1 : MAXindex
        index = MAXindex - i + 1;
        if SORTED_QCLplusTimePoints(index,2) >= CLthrasholdRIGHTtail
            STOPindex = index;
        else 
            break;
        end
    end
    disp(['Number of selected genes in tail is ' num2str(MAXindex-STOPindex)])
    disp(' ')
    for j = STOPindex : STEP : MAXindex
        k = STOPindex + (MAXindex - j) ;
        Q = SORTED_QCLplusTimePoints(k,1);
        CL = SORTED_QCLplusTimePoints(k,2);
        SORTED_Names(k);
        MyTitle = [strjoin(SORTED_Names(k)) ...
            '        Q = ' num2str(Q) ...
            '        CL = ' num2str(CL)];

        disp(['Gene ' num2str(j-STOPindex+1) ' out of ' num2str(MAXindex-STOPindex)])
        if max(SORTED_QCLplusTimePoints(k,1+2:19+2)) >= FilterMaxAmplitude & ...
                mean(SORTED_QCLplusTimePoints(k,1+2:19+2)) >= FilterMean
            AFFIMETRICflag = SORTED_ListAFFIMETRICflags(k,PositionForAFFIflag);
            disp(AFFIMETRICflag)
            if not(FilterAFFIflagsAbsent) | ...
                    not(strcmp(AFFIMETRICflag{1},'YesABSENT'))
                MyFigure = figure;
                subplot(2,1,1)
                plot(Times,SORTED_QCLplusTimePoints(k,1+2:19+2));
                xlabel('time (min)'); ylabel('Gene Expression');
                title(MyTitle);
                subplot(2,1,2)
                autocorr(SORTED_QCLplusTimePoints(k,1+2:19+2))
                waitfor(MyFigure);        
                disp(MyTitle)
                disp(' ')
            end
        else
            disp(['Excluded because of Feng filter'])
            disp(' ')
        end
    end
elseif PassaBanda == 1
    STOPindexMIN = 0;
    for i = 1 : numel(SORTED_QCLplusTimePoints(:,1))
       if SORTED_QCLplusTimePoints(i,2) <= MinCLband
           STOPindexMIN = i;
       else 
           break;
       end
    end
    STOPindexMAX = 0;
    for i = STOPindexMIN : numel(SORTED_QCLplusTimePoints(:,1))
        if SORTED_QCLplusTimePoints(i,2) <= MaxCLband
            STOPindexMAX = i;
        else 
            break;
        end
    end
    disp(['Number of selected genes in tail is ' num2str(STOPindexMAX-STOPindexMIN)])
    disp(' ')
    for k= STOPindexMIN : STEP : STOPindexMAX
        Q = SORTED_QCLplusTimePoints(k,1);
        CL = SORTED_QCLplusTimePoints(k,2);
        SORTED_Names(k);
        MyTitle = [strjoin(SORTED_Names(k)) ...
            '        Q = ' num2str(Q) ...
            '        CL = ' num2str(CL)];

        disp(['Gene ' num2str(k) ' out of ' num2str(STOPindexMAX-STOPindexMIN)])
        if max(SORTED_QCLplusTimePoints(k,1+2:19+2)) >= FilterMaxAmplitude & ...
                mean(SORTED_QCLplusTimePoints(k,1+2:19+2)) >= FilterMean
            AFFIMETRICflag = SORTED_ListAFFIMETRICflags(k,PositionForAFFIflag);
            disp(AFFIMETRICflag)
            if not(FilterAFFIflagsAbsent) | ...
                    not(strcmp(AFFIMETRICflag{1},'YesABSENT'))
                MyFigure = figure;
                subplot(2,1,1)
                plot(Times,SORTED_QCLplusTimePoints(k,1+2:19+2));
                xlabel('time (min)'); ylabel('Gene Expression');
                title(MyTitle);
                subplot(2,1,2)
                autocorr(SORTED_QCLplusTimePoints(k,1+2:19+2))
                waitfor(MyFigure);
                disp(MyTitle)
                disp(' ')
            end
        else
            disp(['Excluded because of Feng filter'])
            disp(' ')
        end
    end
end


