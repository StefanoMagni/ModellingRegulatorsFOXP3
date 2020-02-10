% To write the "probeSetIDs" in a TXT file.

% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 15 May 2017


fileID = fopen('probeSetIDList.txt', 'w');

for i = 1:length(probeSetIDs)
    fprintf(fileID, '%s\n', probeSetIDs{i});
end

fclose(fileID);
