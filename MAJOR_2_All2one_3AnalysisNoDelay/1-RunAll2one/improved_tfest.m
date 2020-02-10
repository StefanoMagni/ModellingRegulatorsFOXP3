function [fitness_tab, models, IO, AIC, gains] = improved_tfest(order, Ts, inputs, outputs)
warning off;
%SystemID using tfest() including offset for SISO systems between lines of data matrix
%Choice between {'iv','svf','gpmf','n4sid','all'}. See doc "tfestoptions"
%opt = tfestOptions('InitMethod','n4sid','Focus',[2*pi/32,2*pi/16]); %Option for system ID ATA initial estimate technique. 
opt = tfestOptions('InitMethod','n4sid'); %Option for system ID ATA initial estimate technique.
%Compute iddata's structure from time series, corresponding model and
%fitness
for i = 1:size(inputs,1)
    in = inputs(i,:);
    parfor j = 1:size(outputs,1) 
        warning off;
        out = outputs(j,:);
        IO{i,j} = iddata(out',[in' ones(size(in'))],Ts); %Include a trivial input for offset
        models{i,j} = tfest(IO{i,j}, [order 0], opt); %0 order model for additional input = offset
        AIC{i,j} = aic(models{i,j});
        [~,fit] = compare(IO{i,j}, models{i,j});
        if fit < 0
            fit = 0;
        end
        g = dcgain(models{i,j}); 
        gains(i,j) = g(1,1);
        if abs(gains(i,j)) <= 0.1 %Test to remove low dcgain models (false systems ID...)
            fitness_tab(i,j) = 0;
        else
            fitness_tab(i,j) = fit;
        end
        %fprintf('Model identified between genes %d and %d, fitness = %d\n',i,j,fit);
    end
    fprintf('%d models identified, %d systems left \n',i*size(inputs,1), size(inputs,1)*size(inputs,1) - i*size(inputs,1));
end

end

