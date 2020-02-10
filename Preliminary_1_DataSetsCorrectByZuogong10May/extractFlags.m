% EXTRACTFLAG - read the CHP data files, sort and assemable into a structure.

% Last modified on 15 May 2017


clear all

% Read CEL raw data and library files
[cels clas exptime] = textread('./newcel.txt','%q%s%s');
numTimePoints = length(clas);   % # of .CEL files; = 74
fp = fopen('probeSetIDList.txt');
probeSetIDList = textscan(fp,'%q');
probeSetIDList = probeSetIDList{1};   % each element is one probeSetID
[numProbeSets, ~] = size(probeSetIDList);

% Path for raw data
celPath = './data/';
libPath = './library/';

% affyread to load raw data
% % TimeSeriesMatrix = zeros(numProbeSets, numTimePoints);
FlagsCellArray = cell(numProbeSets, numTimePoints);
logFileID = fopen('process.log', 'w');
for i = 1:numTimePoints
    % celStruct = affyread('./data/GSM723834.CEL');
    % maimage(celStruct, 'Intensity')
    % probeStruct = celintensityread(cels,'HG-U133_Plus_2.cdf', 'celpath', celPath, 'cdfpath', libPath);
    chpStruct_ti = affyread([celPath cels{i} '.CHP'], libPath);
    if chpStruct_ti.NumProbeSets ~= numProbeSets
        % To check that this CHP file contains the info of all genes (probeSets)
        disp('Warning: records of some probe sets are missing!');
    end
    fprintf(logFileID, '>> Data File: %s (%i/%i):\n', cels{i}, i, numTimePoints);
    
    for j = 1:numProbeSets
        probeSetID = probeSetIDList{j};
        PPStruct_ti = probesetlookup(chpStruct_ti, probeSetID);
        geneIndex = PPStruct_ti.CDFIndex;
% %         TimeSeriesMatrix(j, i) = chpStruct_ti.ProbeSets(geneIndex).Signal;
        FlagsCellArray{j,i} = chpStruct_ti.ProbeSets(geneIndex).Detection;
        % {j,i}: the j-th probeSet (gene) of the i-th time point (.CEL, .CHP)
        
        if rem(j, 500) == 0
            fprintf(logFileID, '#ProbeSet Processed: %i/%i\n', j, numProbeSets);
        end
    end
end
fclose(logFileID);

% Build the original record in the DataMatrix type
% % timeIndexInfo = strcat(clas, ',', exptime);
% % import bioma.data.*
% % TimeSeriesRawNomalization = ...
% %     DataMatrix(TimeSeriesMatrix, probeSetIDList, timeIndexInfo);
% % % TimeSeriesrawnomalization: the time series used in Feng's PLAU paper


% Save variables
savePath = 'correctedDataMat/';
save([savePath './FlagsCellArray.mat'], 'FlagsCellArray');
% % save([savePath './TimeSeriesRawNomalization.mat'], 'TimeSeriesRawNomalization');


% End 15 May 2017
