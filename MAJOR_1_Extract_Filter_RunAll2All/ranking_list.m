function ranks = ranking_list(fit_table,names,AIC)
%Return table and names of inputs/outputs ranks from fitness_t (i fit table
%= inputs; y fit table = outputs)
dim1 = size(fit_table,1);
B = reshape(fit_table,[1 dim1*dim1]);
[~,idx] = sort(B,'descend'); 

for i = 1:(dim1*dim1)
    ranks{i,1} = fit_table(idx(i)); %Fitness value
    %Input
    if mod(idx(i),dim1) == 0
        ranks{i,2} = names(dim1);
    else
        ranks{i,2} = names(mod(idx(i),dim1)); 
    end
    ranks{i,3} = names(ceil(idx(i)/dim1)); %Output
    ranks{i,4} = AIC(idx(i)); % Corresponding AIC Value
end

end