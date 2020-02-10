%% To run all2one
% Main script written by Laurent on Dec, 2017
% Modified by Rucha S. on Dec, 2017

function [fitness_tab, models, IO, AIK, dcs] = all_to_one_delay_tfest_MOD11Dec(inputs,output,order,delay)
%SystemID using tfest() including offset for SISO systems between lines of data matrix
init_opt = 'n4sid'; %Choice between {'iv','svf','gpmf','n4sid','all'}. See doc "tfestoptions"
%opt = tfestOptions('InitMethod','n4sid','Focus',[2*pi/32,2*pi/16]); %Option for system ID ATA initial estimate technique. 
opt = tfestOptions('InitMethod','n4sid'); %Option for system ID ATA initial estimate technique.
warning off;
%Compute iddata's structure from time series, corresponding model and
%fitness
%Compute iddata's structure from time series, corresponding model and
%fitness
models = cell(1,size(inputs,1));
fitness_tab = zeros(1,size(inputs,1));
AIK = zeros(1,size(inputs,1));
dcs = zeros(1,size(inputs,1));
errorA = NaN(1,size(inputs,1));

parfor i = 1:size(inputs,1)
    warning off;    
    %Model computation
    IO{i} = iddata(output',[inputs(i,:)' ones(size(output'))],1); %Include a trivial input for offset
    
    try
        models{i} = tfest(IO{i},[order 0], opt, 'InputDelay', [delay 0]); %0 order model for additional input = offset
        errorA(i) = 0;
    catch
        errorA(i) = 1;
    end
    if errorA(i) == 0
        [~,fit] = compare(IO{i}, models{i});
        AIK(i) = aic(models{i},'AICc');

        %Input sufficiency check
        gain = dcgain(models{i}); 
        g = gain(1,1);
        dcs(i) = g;
        if order == 1
            if abs(g) <= 0.1 
                fitness_tab(i) = 0;
            else
                fitness_tab(i) = fit;
            end 
        else
            fitness_tab(i) = fit;
        end
        if fitness_tab(i) < 0
            fitness_tab(i) = 0;
        end
    else
        fitness_tab(i) = NaN;
        dcs(i) = NaN;
        AIK(i) = NaN;
        models{i} = NaN;
    end

    fprintf('Model from gene %d computed\n',i);
end

end

