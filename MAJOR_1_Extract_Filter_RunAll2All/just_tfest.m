function [fitness_tab, models, IO, AIK, gain] = just_tfest(order, Ts, data, opt)
%SystemID using tfest() including offset for SISO systems between lines of data matrix
warning off;
%Compute iddata's structure from time series, corresponding model and
%fitness
for i = 1:size(data,1)
    parfor j = 1:size(data,1)
        warning off;
        if i ~= j
            IO{i,j} = iddata(data(j,:)',[data(i,:)' ones(size(data(i,:)'))],Ts); %Include a trivial input for offset
            models{i,j} = tfest(IO{i,j}, [order 0], opt); %0 order model for additional input = offset
            [~,fit] = compare(IO{i,j}, models{i,j});
            AIK(i,j) = aic(models{i,j});
            gain{i,j} = dcgain(models{i,j}); 
%             if abs(gain) <= 0.1 %Test to remove low dcgain models (false systems ID...)
%                 fitness_tab(i,j) = 0;
%             else
            fitness_tab(i,j) = fit;
            if fitness_tab(i,j) < 0
                fitness_tab(i,j) = 0;
            end
%             end
            %fprintf('Model identified between genes %d and %d, fitness = %d\n',i,j,fit);
        else
            IO{i,j} = NaN;
            models{i,j} = NaN;
            fitness_tab(i,j) = 0;
            AIK(i,j) = 0;
            gain{i,j} = 0;
        end
    end
    fprintf('%d models identified, %d systems left \n',i*size(data,1), size(data,1)*size(data,1) - i*size(data,1));
end

end

