
function [ranks, idx] = ranking_list_FOXP3(fit_score,names,ActRep)
%Return table and names of inputs/outputs ranks from fitness_t (i fit table
%= inputs; y fit table = outputs)
[~,idx] = sort(fit_score,'descend'); 

for i = 1:size(fit_score,1)
    ranks{i,1} = idx(i);
    ranks{i,2} = fit_score(idx(i)); % Fitness value
    ranks{i,3} = names{idx(i)};     % Gene Names
    ranks{i,4} = ActRep(idx(i));    % Activation-Repression
end

end












% function ranks = ranking_list_FOXP3(fit_table,names,AIC,dcs,Loc)
% %Return table and names of inputs/outputs ranks from fitness_t (i fit table
% %= inputs; y fit table = outputs)
% [~,idx] = sort(fit_table,'descend'); 
% 
% for i = 1:size(fit_table,1)
%     ranks{i,1} = idx(i);
%     ranks{i,2} = Loc(i);
%     ranks{i,3} = fit_table(idx(i)); %Fitness value
%     ranks{i,4} = names{idx(i)}; %Input
%     ranks{i,5} = dcs(idx(i)); %dcgqin
%     ranks{i,6} = AIC(idx(i)); % Corresponding AIC Value
% end
% 
% end