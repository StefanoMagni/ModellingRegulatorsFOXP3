% Sript from LAURENT MOMBAERTS to plot a graph of the results from All2All
% algorithm. Modified by STEFANO MAGNI

%========================================================================
%                             OLD 
%========================================================================


function plot_biographMODforVectors(data, threshold, species, FigName, dcgains, mRNA_name_SingleElement)%,dcgains)

%             data=fitness;
%             species=mRNAs_names;
            
    species(1) = {mRNA_name_SingleElement};

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
    cmMOD = [cm,zeros(size(cm,1),size(cm,1)-1)];
    dcgainsMOD = [dcgains,zeros(size(dcgains,1),size(dcgains,1)-1)];
    %cmMOD = [cm;zeros(size(cm,2)-1,size
    ids = species;
    bg2 = biograph(cmMOD,ids,'ShowWeights', 'on');%'ArrowSize', 5, 'EdgeFontSize', 15);
    %bg2 = biograph(cm,ids);%,
    get(bg2.nodes,'ID');
    %view(bg2);

    %gene_name    = cell(1, 41);
    %gene_name(:) = {'FOXP3'};
    for i = 1:size(data,1)
        for j = 1:size(data,2)%,2)
            if cmMOD(i,j) > 0
                %If a link, check which type of regulation
                if dcgainsMOD(i,j) >= 0
                    %Blue
                    %size(species(j));
                    %size(mRNA_name_SingleElement);

                    %disp(species(j));
                    set(getedgesbynodeid(bg2,species(i),species(j)),'LineColor',[0 0 1]);
                    %set(getedgesbynodeid(bg2,species(i),species(j)),'LineColor',[0 0 1]);
                else
                    %red

                    %disp(species(j));
                    set(getedgesbynodeid(bg2,species(i),species(j)),'LineColor',[1 0 0]);
                    %set(getedgesbynodeid(bg2,species(i),species(j)),'LineColor',[1 0 0]);
                end
            end
        end
    end

    %view(bg2);

    g = biograph.bggui(bg2);
    f = figure;
    copyobj(g.biograph.hgAxes,f);
    %savefig(FigName);
    temp=[FigName,'.eps']; saveas(gcf,temp,'epsc');

end







%========================================================================
%                          NEW 
%========================================================================

% % 
% % function plot_biographMODforVectors(data, threshold, species, FigName, ...
% %     dcgains, mRNA_name_SingleElement)%,dcgains)
% % 
% % 	data = [0;data];
% % 	species = [mRNA_name_SingleElement ; species];
% %     dcgains = [0 ; dcgains];
% %             
% % %     species(1) = {mRNA_name_SingleElement};
% % 
% %     for i = 1:size(data,1) 
% %         for j = 1:size(data,2) 
% %             if data(i,j) >= threshold
% %                 net(i,j) = floor(data(i,j)); 
% %             else
% %                 net(i,j) = 0; 
% %             end
% %         end
% %     end
% % 
% %     cm = net;
% %     cmMOD = [cm,zeros(size(cm,1),size(cm,1)-1)];
% %     dcgainsMOD = [dcgains,zeros(size(dcgains,1),size(dcgains,1)-1)];
% %     %cmMOD = [cm;zeros(size(cm,2)-1,size
% %     ids = species;
% %     bg2 = biograph(cmMOD,ids,'ShowWeights', 'on');%'ArrowSize', 5, 'EdgeFontSize', 15);
% %     %bg2 = biograph(cm,ids);%,
% %     get(bg2.nodes,'ID');
% %     %view(bg2);
% % 
% % %     % If a Gene Name appear twice add probe ID
% % %     [un idx_last idx] = unique(species(:,1));
% % %     unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {sort(x)});
% % %     
% % %     for k = 1:size(unique_idx,1)
% % %         if size(unique_idx{k},1) >= 2
% % %            speciesNew = [sprintf(species{unique_idx{k}}),'-',num2str(k)]; 
% % %            species{k+1} = speciesNew;
% % %         end
% % %     end
% % 
% %     %cmMOD(1,1)=0;
% %     %gene_name    = cell(1, 41);
% %     %gene_name(:) = {'FOXP3'};
% %     for i = 1:size(data,1)
% %         for j = 1:size(data,2)%,2)
% %             if cmMOD(i,j) > 0
% %                 %If a link, check which type of regulation
% %                 if dcgainsMOD(i,j) >= 0
% %                     %Blue
% %                     %size(species(j));
% %                     %size(mRNA_name_SingleElement);
% % 
% %                     disp(i);
% %                     disp(j);
% %                     disp(species(j));
% % 
% %                     set(getedgesbynodeid(bg2,species(i),species(j)),'LineColor',[0 0 1]);
% %                     %set(getedgesbynodeid(bg2,species(i),species(j)),'LineColor',[0 0 1]);
% %                 else
% %                     %red
% % 
% %                     %disp(species(j));
% %                     set(getedgesbynodeid(bg2,species(i),species(j)),'LineColor',[1 0 0]);
% %                     %set(getedgesbynodeid(bg2,species(i),species(j)),'LineColor',[1 0 0]);
% %                 end
% %             end
% %         end
% %     end
% % 
% %     %view(bg2);
% % 
% %     g = biograph.bggui(bg2);
% %     f = figure;
% %     copyobj(g.biograph.hgAxes,f);
% %     %savefig(FigName);
% %     temp=[FigName,'.eps']; saveas(gcf,temp,'epsc');
% % 
% % end

