%% To plot the PLAU inferred.. Feng et al. 2012, Figure 2 - Teff and Treg cell data - but in log scale.
% Version: 19/01/2017 Rucha S.

%% COMPARE AND ALIGN PROBE ID+OTHER INFO W.R.T. THAT IN Known_resp_genes.mat

clc;
clear;
close all;

% load('Known_resp_genes.mat')     % Note: Row 1 of title is removed from coresponding Excel file as well.
% load('expr_info_data.mat')       % Note: Row 1 of title is removed from coresponding Excel file as well.
% load('expr_log_info.mat')   

load('known_PLAU_log.mat')

%------------------
    a = 'IL2';
    b = 'IL4';
    c = 'IL5';
    d = 'CSF2';
    e = 'IL13';
    f = 'IFNG';
    g = 'CCL20';
    h = 'FOXP3';
    k = 'LRRC32';
    
    match_set = {a,b,c,d,e,f,g,h,k}';  % Set of genes to find & plot

%-----------------------
%%Uncomment this part in case one wants to search from the big, main file!
% count = 0;
% for i=1:9
%     count = count+1;    
%     % 'find_exact' search for each Gene Name of 'Known_resp_genes' in 'expr_info_data', 
%     % and provides location of that exact matching string in 'expr_info_data'
%     
%     strs_PLAU_log = match_set(count,1);
%     find_PLAU_log = find(ismember(expr_log_info(1:54675,4),strs_PLAU_log));
%     loc_PLAU_log(count,1:length(find_PLAU_log)) = find_PLAU_log; 
%     
%     loc_PLAU_log_row = loc_PLAU_log';
%     non_zero_PLAU_log = loc_PLAU_log_row(loc_PLAU_log_row~=0);
%     
%     
% end  
% 
% known_PLAU_log = expr_log_info(non_zero_PLAU_log(:,:),:);
%-----------------------


PLAU_log_double = cell2mat(known_PLAU_log(:,6:81));



%% PLOTS

T = 360;
t_total = 0:20:T;


figure(1)
grid on;
title('Known Responsive Genes') 

subplot(3,2,1)
plot(t_total,PLAU_log_double(1,39:57),'bd-','MarkerFaceColor','b','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(1,58:76),'ms-','MarkerFaceColor','m','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(1,1:19),'y^-','LineWidth',1.5);   hold on;
plot(t_total,PLAU_log_double(1,20:38),'cx-','LineWidth',1.5);  hold on;
title('IL2')
xlabel('Time (min)'); ylabel('mRNA exprn. level')
legend('Treg-1','Treg-2','Teff-1','Teff-2') 
% axis([0 360 0 20000])

subplot(3,2,2)
plot(t_total,PLAU_log_double(2:3,39:57),'bd-','MarkerFaceColor','b','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(2:3,58:76),'ms-','MarkerFaceColor','m','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(2:3,1:19),'y^-','LineWidth',1.5);   hold on;
plot(t_total,PLAU_log_double(2:3,20:38),'cx-','LineWidth',1.5);  hold on;
title('IL4')
xlabel('Time (min)'); ylabel('mRNA exprn. level')
legend('Treg-1','Treg-2','Teff-1','Teff-2') 
% axis([0 360 0 30000])

subplot(3,2,3)
plot(t_total,PLAU_log_double(4,39:57),'bd-','MarkerFaceColor','b','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(4,58:76),'ms-','MarkerFaceColor','m','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(4,1:19),'y^-','LineWidth',1.5);   hold on;
plot(t_total,PLAU_log_double(4,20:38),'cx-','LineWidth',1.5);  hold on;
title('IL5')
xlabel('Time (min)'); ylabel('mRNA exprn. level')
legend('Treg-1','Treg-2','Teff-1','Teff-2') 
% axis([0 360 0 25000])

subplot(3,2,4) % row 5 of CSF2 not included in the plot
plot(t_total,PLAU_log_double(6,39:57),'bd-','MarkerFaceColor','b','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(6,58:76),'ms-','MarkerFaceColor','m','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(6,1:19),'y^-','LineWidth',1.5);   hold on;
plot(t_total,PLAU_log_double(6,20:38),'cx-','LineWidth',1.5);  hold on;
title('CSF2')
xlabel('Time (min)'); ylabel('mRNA exprn. level')
legend('Treg-1','Treg-2','Teff-1','Teff-2') 
% axis([0 360 0 40000])

subplot(3,2,5)
plot(t_total,PLAU_log_double(7,39:57),'bd-','MarkerFaceColor','b','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(7,58:76),'ms-','MarkerFaceColor','m','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(7,1:19),'y^-','LineWidth',1.5);   hold on;
plot(t_total,PLAU_log_double(7,20:38),'cx-','LineWidth',1.5);  hold on;
title('IL13')
xlabel('Time (min)'); ylabel('mRNA exprn. level')
legend('Treg-1','Treg-2','Teff-1','Teff-2') 
% axis([0 360 0 40000])

subplot(3,2,6)
plot(t_total,PLAU_log_double(8,39:57),'bd-','MarkerFaceColor','b','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(8,58:76),'ms-','MarkerFaceColor','m','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(8,1:19),'y^-','LineWidth',1.5);   hold on;
plot(t_total,PLAU_log_double(8,20:38),'cx-','LineWidth',1.5);  hold on;
title('IFNG')
xlabel('Time (min)'); ylabel('mRNA exprn. level')
legend('Treg-1','Treg-2','Teff-1','Teff-2') 
% axis([0 360 0 15000])


%--------------------------




figure(2)
grid on;
title('Known Responsive Genes') 

subplot(1,3,1)
plot(t_total,PLAU_log_double(9,39:57),'bd-','MarkerFaceColor','b','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(9,58:76),'ms-','MarkerFaceColor','m','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(9,1:19),'y^-','LineWidth',1.5);   hold on;
plot(t_total,PLAU_log_double(9,20:38),'cx-','LineWidth',1.5);  hold on;
title('CCL20')
xlabel('Time (min)'); ylabel('mRNA exprn. level')
legend('Treg-1','Treg-2','Teff-1','Teff-2') 
% axis([0 360 0 6000])

subplot(1,3,2) % row 10-11 for FOXP3 not included in the plot
plot(t_total,PLAU_log_double(12,39:57),'bd-','MarkerFaceColor','b','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(12,58:76),'ms-','MarkerFaceColor','m','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(12,1:19),'y^-','LineWidth',1.5);   hold on;
plot(t_total,PLAU_log_double(12,20:38),'cx-','LineWidth',1.5);  hold on;
title('FOXP3')
xlabel('Time (min)'); ylabel('mRNA exprn. level')
legend('Treg-1','Treg-2','Teff-1','Teff-2') 
% axis([0 360 0 1000])

subplot(1,3,3)
plot(t_total,PLAU_log_double(13,39:57),'bd-','MarkerFaceColor','b','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(13,58:76),'ms-','MarkerFaceColor','m','LineWidth',1.5);  hold on;
plot(t_total,PLAU_log_double(13,1:19),'y^-','LineWidth',1.5);   hold on;
plot(t_total,PLAU_log_double(13,20:38),'cx-','LineWidth',1.5);  hold on;
title('LRRC32')
xlabel('Time (min)'); ylabel('mRNA exprn. level')
legend('Treg-1','Treg-2','Teff-1','Teff-2') 
% axis([0 360 0 6000])





%     *****************       END       *****************




    
    
    
    

