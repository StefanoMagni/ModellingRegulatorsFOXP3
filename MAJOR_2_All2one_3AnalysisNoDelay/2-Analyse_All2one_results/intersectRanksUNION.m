%% Function


function [Intersection, ProbeNumberIntersection] = intersectRanks(ranks1,ranks2, ...
    ColumnNumberForIntersect, ProbeNumber)

    N = ColumnNumberForIntersect;

    for ii = 1:size(ranks1,1)
        loc(ii,:) = strcmp(ranks1(ii,N),ranks2(:,N));
    end

    [r c] = find(loc == 1);
    Intersection = [ranks1(r,:) , ranks2(c,:)];
    ProbeNumberIntersection = [ProbeNumber(r,1),ProbeNumber(c,2)];
    
end



