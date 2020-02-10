eclear

% In the "newcel.txt", we remove "GSM723850" (Teff-1 320min)
% Note that Treg-1 120min is missed. We need to interpolate two points.
[cels clas exptime] = textread('newcel.txt','%q%u%s');
n=length(clas);
celPath = './data';
libPath = './library';

%% Normalization from Raw Data (*.CEL) 
% GCRMA
expr_gcrma = affygcrma(cels,'HG-U133_Plus_2.cdf', 'HG-U133_Plus_2.probe_tab',...
    'celpath', celPath, 'cdfpath', libPath, 'seqpath', libPath);
eObjDataMatrix = 2.^expr_gcrma;

% % GCRMA in linear
% eObjDataMatrix = affygcrma(cels,'HG-U133_Plus_2.cdf', 'HG-U133_Plus_2.probe_tab',...
%     'celpath', celPath, 'cdfpath', libPath, 'seqpath', libPath, 'Output', ...
%                        'linear');

%% Save Nomalized Results
%   expr_gcrma: the normalized data by gcrma in log2 scale (type: DataMatrix)
%   e: .... in linear scale (type: double)
e = double(eObjDataMatrix);
savePath = './correctedMatData/';
save([savePath 'expr_gcrma.mat'], 'expr_gcrma');
save([savePath 'e.mat'], 'e');

%% Generate "en.mat", which is the "e.mat" after linear interpolation
% You may use other interpolation here.
en = [e(:, 1:16) ...
      (e(:, 16)+e(:, 17))/2  ...
      e(:, 17:43) ...
      (e(:, 43) + e(:, 44))/2 ...
      e(:, 44:74)];
probeSetIDs = expr_gcrma.RowNames;
save([savePath 'en.mat'], 'en', 'probeSetIDs');


% Notes of DATA:
% + data[1:75]: teff1, teff2, treg1, treg2; treg1 120min missing.
% + annotation: "HG-U133_Plus_2.annot.revised.csv"; 