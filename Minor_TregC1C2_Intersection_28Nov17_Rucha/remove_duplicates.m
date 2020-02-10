% TO REMOVE ANY GENES WITH DUPLICATE LOCATION AND GENE NAME


function [withoutduplicates] = remove_duplicates(withduplicates)

    wd=withduplicates;
    [~,idx]=unique( strcat(wd(:,1),wd(:,3)) );
    withoutduplicates=wd(idx,:);



end



