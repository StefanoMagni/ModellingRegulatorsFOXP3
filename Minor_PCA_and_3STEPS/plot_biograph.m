% Sript from LAURENT MOMBAERTS to plot a graph of the results from All2All
% algorithm. Modified by STEFANO MAGNI

function plot_biograph(data,threshold,species, FigName, dcgains)%,dcgains)

for i = 1:size(data,1) 
    for j = 1:size(data,2) 
        if data(i,j) >= threshold
            net(i,j) = floor(data(i,j)); 
        else
            net(i,j) = 0; 
        end
    end
end

cm = net;
ids = species;
bg2 = biograph(cm,ids,'ShowWeights', 'on');%'ArrowSize', 5, 'EdgeFontSize', 15);
%bg2 = biograph(cm,ids);%,
get(bg2.nodes,'ID');
%view(bg2);

for i = 1:size(data,1)
    for j = 1:size(data,1)%,2)
        if net(i,j) > 0
            %If a link, check which type of regulation
            if dcgains{i,j} >= 0
                %Blue
                set(getedgesbynodeid(bg2,species(i),species(j)),'LineColor',[0 0 1]);
            else
                %red
                set(getedgesbynodeid(bg2,species(i),species(j)),'LineColor',[1 0 0]);
            end
        end
    end
end

%view(bg2);

g = biograph.bggui(bg2);
f = figure;
copyobj(g.biograph.hgAxes,f);
savefig(FigName);

end