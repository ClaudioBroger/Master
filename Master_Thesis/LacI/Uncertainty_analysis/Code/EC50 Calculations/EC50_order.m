%% Calculate order of EC50 and differences in parameters 

%%Hyperspace

load('04-Apr-2023ec50_1_random_200.mat');
load('04-Apr-2023ec50_2_random_200.mat');
load('04-Apr-2023ec50_3_random_200.mat');

ec50_1over2 = find(ec50_1(:,2) >= ec50_2(:,2));
ec50_2over3 = find(ec50_2(:,2) >= ec50_3(:,2));
ec50_1over3 = find(ec50_1(:,2) >= ec50_3(:,2));

%No EC50 in the wrong order found for hyperspace

%%PESTO
load('03-Apr-2023ec50_1_random_pesto_200.mat');
load('03-Apr-2023ec50_2_random_pesto_200.mat');
load('03-Apr-2023ec50_3_random_pesto_200.mat');

ec50_1over2_pesto = find(ec50_1_pesto(:,2) >= ec50_2_pesto(:,2));
ec50_2over3_pesto = find(ec50_2_pesto(:,2) >= ec50_3_pesto(:,2));
ec50_1over3_pesto = find(ec50_1_pesto(:,2) >= ec50_3_pesto(:,2));

%15 EC50 from first repression coefficient found to be over EC50 from
%second repression coefficient

load('03-Apr-2023_rand_parameter_Pesto.mat');

rand_parameter_wrong_order = rand_parameter(ec50_1over2_pesto,:);

rand_parameter_right_order = rand_parameter;
rand_parameter_right_order(ec50_1over2_pesto,:) = [];

rand_parameter_right_order = rand_parameter_right_order(randi([1 height(rand_parameter_right_order)],1,height(rand_parameter_wrong_order)),:);

%ANOVA for each parameter
rand_parameter_right_order.dCit = [];
rand_parameter_right_order.mu = [];
rand_parameter_right_order.kmaturation = [];
rand_parameter_right_order.nMperUnit = [];
rand_parameter_wrong_order.dCit = [];
rand_parameter_wrong_order.mu = [];
rand_parameter_wrong_order.kmaturation = [];
rand_parameter_wrong_order.nMperUnit = [];

kLacI_pesto = [table2array(rand_parameter_right_order(:,"kLacI")) table2array(rand_parameter_wrong_order(:,"kLacI"))];
group = ["right" "wrong"];
[p_kLacI_pesto, tbl_kLacI_pesto, stats_kLacI_pesto] = anova1(kLacI_pesto, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO - kLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_kLacI_pesto);
hold on
title('Multcompare order of EC50 PESTO - kLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_kLacI_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_kLacI_PESTO', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_kLacI_PESTO', '.jpg']);

close all 


kCit_pesto = [table2array(rand_parameter_right_order(:,"kCit")) table2array(rand_parameter_wrong_order(:,"kCit"))];
group = ["right" "wrong"];
[p_kCit_pesto, tbl_kCit_pesto, stats_kCit_pesto] = anova1(kCit_pesto, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO - kCit', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_kCit_pesto);
hold on
title('Multcompare order of EC50 PESTO - kCit', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_kCit_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_kCit_PESTO', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_kCit_PESTO', '.jpg']);
close all  

dLacI_pesto = [table2array(rand_parameter_right_order(:,"dLacI")) table2array(rand_parameter_wrong_order(:,"dLacI"))];
group = ["right" "wrong"];
[p_dLacI_pesto, tbl_dLacI_pesto, stats_dLacI_pesto] = anova1(dLacI_pesto, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO - dLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_dLacI_pesto);
hold on
title('Multcompare order of EC50 PESTO - dLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_dLacI_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_dLacI_PESTO', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_dLacI_PESTO', '.jpg']);

close all  

LacIrep_pesto = [table2array(rand_parameter_right_order(:,"LacIrep")) table2array(rand_parameter_wrong_order(:,"LacIrep"))];
group = ["right" "wrong"];
[p_LacIrep_pesto, tbl_LacIrep_pesto, stats_LacIrep_pesto] = anova1(LacIrep_pesto, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO - LacIrep', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_LacIrep_pesto);
hold on
title('Multcompare order of EC50 PESTO - LacIrep', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep_PESTO', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_LacIrep_PESTO', '.jpg']);

close all  

nLacI_pesto = [table2array(rand_parameter_right_order(:,"nLacI")) table2array(rand_parameter_wrong_order(:,"nLacI"))];
group = ["right" "wrong"];
[p_nLacI_pesto, tbl_nLacI_pesto, stats_nLacI_pesto] = anova1(nLacI_pesto, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO - nLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_nLacI_pesto);
hold on
title('Multcompare order of EC50 PESTO - nLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_nLacI_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_nLacI_PESTO', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_nLacI_PESTO', '.jpg']);

close all  

CitL_pesto = [table2array(rand_parameter_right_order(:,"CitL")) table2array(rand_parameter_wrong_order(:,"CitL"))];
group = ["right" "wrong"];
[p_CitL_pesto, tbl_CitL_pesto, stats_CitL_pesto] = anova1(CitL_pesto, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO - CitL', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_CitL_pesto);
hold on
title('Multcompare order of EC50 PESTO - CitL', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_CitL_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_CitL_PESTO', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_CitL_PESTO', '.jpg']);

close all  

KdLacI_pesto = [table2array(rand_parameter_right_order(:,"KdLacI")) table2array(rand_parameter_wrong_order(:,"KdLacI"))];
group = ["right" "wrong"];
[p_KdLacI_pesto, tbl_KdLacI_pesto, stats_KdLacI_pesto] = anova1(KdLacI_pesto, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO - KdLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_KdLacI_pesto);
hold on
title('Multcompare order of EC50 PESTO - KdLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_KdLacI_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_KdLacI_PESTO', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_KdLacI_PESTO', '.jpg']);

close all 

LacIrep2_pesto = [table2array(rand_parameter_right_order(:,"LacIrep2")) table2array(rand_parameter_wrong_order(:,"LacIrep2"))];
group = ["right" "wrong"];
[p_LacIrep2_pesto, tbl_LacIrep2_pesto, stats_LacIrep2_pesto] = anova1(LacIrep2_pesto, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO - LacIrep2', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_LacIrep2_pesto);
hold on
title('Multcompare order of EC50 PESTO - LacIrep2', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep2_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep2_PESTO', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_LacIrep2_PESTO', '.jpg']);

close all 

LacIrep3_pesto = [table2array(rand_parameter_right_order(:,"LacIrep3")) table2array(rand_parameter_wrong_order(:,"LacIrep3"))];
group = ["right" "wrong"];
[p_LacIrep3_pesto, tbl_LacIrep3_pesto, stats_LacIrep3_pesto] = anova1(LacIrep3_pesto, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO - LacIrep3', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_LacIrep3_pesto);
hold on
title('Multcompare order of EC50 PESTO - LacIrep3', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep3_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep3_PESTO', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_LacIrep3_PESTO', '.jpg']);

close all 


%%PESTO CI
load('04-Apr-2023ec50_1_random_pesto_CI_200.mat');
load('04-Apr-2023ec50_2_random_pesto_CI_200.mat');
load('04-Apr-2023ec50_3_random_pesto_CI_200.mat');

ec50_1over2_pesto = find(ec50_1_pesto(:,2) >= ec50_2_pesto(:,2));
ec50_2over3_pesto = find(ec50_2_pesto(:,2) >= ec50_3_pesto(:,2));
ec50_1over3_pesto = find(ec50_1_pesto(:,2) >= ec50_3_pesto(:,2));

%15 EC50 from first repression coefficient found to be over EC50 from
%second repression coefficient

load('03-Apr-2023_rand_parameter_Pesto_cI-200.mat');

rand_parameter_wrong_order = rand_parameter(ec50_1over2_pesto,:);

rand_parameter_right_order = rand_parameter;
rand_parameter_right_order(ec50_1over2_pesto,:) = [];

rand_parameter_right_order = rand_parameter_right_order(randi([1 height(rand_parameter_right_order)],1,height(rand_parameter_wrong_order)),:);

%ANOVA for each parameter
rand_parameter_right_order.dCit = [];
rand_parameter_right_order.mu = [];
rand_parameter_right_order.kmaturation = [];
rand_parameter_right_order.nMperUnit = [];
rand_parameter_wrong_order.dCit = [];
rand_parameter_wrong_order.mu = [];
rand_parameter_wrong_order.kmaturation = [];
rand_parameter_wrong_order.nMperUnit = [];

kLacI_pesto_CI = [table2array(rand_parameter_right_order(:,"kLacI")) table2array(rand_parameter_wrong_order(:,"kLacI"))];
group = ["right" "wrong"];
[p_kLacI_pesto_CI, tbl_kLacI_pesto_CI, stats_kLacI_pesto_CI] = anova1(kLacI_pesto_CI, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO CI - kLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_kLacI_pesto_CI);
hold on
title('Multcompare order of EC50 PESTO CI - kLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_kLacI_table_PESTO_CI', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_kLacI_PESTO_CI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_kLacI_PESTO_CI', '.jpg']);

close all 


kCit_pesto_CI = [table2array(rand_parameter_right_order(:,"kCit")) table2array(rand_parameter_wrong_order(:,"kCit"))];
group = ["right" "wrong"];
[p_kCit_pesto_CI, tbl_kCit_pesto_CI, stats_kCit_pesto_CI] = anova1(kCit_pesto_CI, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO CI - kCit', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_kCit_pesto_CI);
hold on
title('Multcompare order of EC50 PESTO CI - kCit', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_kCit_table_PESTO_CI', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_kCit_PESTO_CI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_kCit_PESTO_CI', '.jpg']);
close all  

dLacI_pesto_CI = [table2array(rand_parameter_right_order(:,"dLacI")) table2array(rand_parameter_wrong_order(:,"dLacI"))];
group = ["right" "wrong"];
[p_dLacI_pesto_CI, tbl_dLacI_pesto_CI, stats_dLacI_pesto_CI] = anova1(dLacI_pesto_CI, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO CI - dLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_dLacI_pesto_CI);
hold on
title('Multcompare order of EC50 PESTO CI - dLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_dLacI_table_PESTO_CI', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_dLacI_PESTO_CI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_dLacI_PESTO_CI', '.jpg']);

close all  

LacIrep_pesto_CI = [table2array(rand_parameter_right_order(:,"LacIrep")) table2array(rand_parameter_wrong_order(:,"LacIrep"))];
group = ["right" "wrong"];
[p_LacIrep_pesto_CI, tbl_LacIrep_pesto_CI, stats_LacIrep_pesto_CI] = anova1(LacIrep_pesto_CI, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO CI - LacIrep', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_LacIrep_pesto_CI);
hold on
title('Multcompare order of EC50 PESTO CI - LacIrep', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep_table_PESTO_CI', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep_PESTO_CI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_LacIrep_PESTO_CI', '.jpg']);

close all  

nLacI_pesto_CI = [table2array(rand_parameter_right_order(:,"nLacI")) table2array(rand_parameter_wrong_order(:,"nLacI"))];
group = ["right" "wrong"];
[p_nLacI_pesto_CI, tbl_nLacI_pesto_CI, stats_nLacI_pesto_CI] = anova1(nLacI_pesto_CI, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO CI - nLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_nLacI_pesto_CI);
hold on
title('Multcompare order of EC50 PESTO CI - nLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_nLacI_table_PESTO_CI', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_nLacI_PESTO_CI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_nLacI_PESTO_CI', '.jpg']);

close all  

CitL_pesto_CI = [table2array(rand_parameter_right_order(:,"CitL")) table2array(rand_parameter_wrong_order(:,"CitL"))];
group = ["right" "wrong"];
[p_CitL_pesto_CI, tbl_CitL_pesto_CI, stats_CitL_pesto_CI] = anova1(CitL_pesto_CI, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO CI - CitL', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_CitL_pesto_CI);
hold on
title('Multcompare order of EC50 PESTO CI - CitL', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_CitL_table_PESTO_CI', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_CitL_PESTO_CI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_CitL_PESTO_CI', '.jpg']);

close all  

KdLacI_pesto_CI = [table2array(rand_parameter_right_order(:,"KdLacI")) table2array(rand_parameter_wrong_order(:,"KdLacI"))];
group = ["right" "wrong"];
[p_KdLacI_pesto_CI, tbl_KdLacI_pesto_CI, stats_KdLacI_pesto_CI] = anova1(KdLacI_pesto_CI, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO CI - KdLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_KdLacI_pesto_CI);
hold on
title('Multcompare order of EC50 PESTO CI - KdLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_KdLacI_table_PESTO_CI', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_KdLacI_PESTO_CI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_KdLacI_PESTO_CI', '.jpg']);

close all 

LacIrep2_pesto_CI = [table2array(rand_parameter_right_order(:,"LacIrep2")) table2array(rand_parameter_wrong_order(:,"LacIrep2"))];
group = ["right" "wrong"];
[p_LacIrep2_pesto_CI, tbl_LacIrep2_pesto_CI, stats_LacIrep2_pesto_CI] = anova1(LacIrep2_pesto_CI, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO CI - LacIrep2', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_LacIrep2_pesto_CI);
hold on
title('Multcompare order of EC50 PESTO CI - LacIrep2', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep2_table_PESTO_CI', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep2_PESTO_CI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_LacIrep2_PESTO_CI', '.jpg']);

close all 

LacIrep3_pesto_CI = [table2array(rand_parameter_right_order(:,"LacIrep3")) table2array(rand_parameter_wrong_order(:,"LacIrep3"))];
group = ["right" "wrong"];
[p_LacIrep3_pesto_CI, tbl_LacIrep3_pesto_CI, stats_LacIrep3_pesto_CI] = anova1(LacIrep3_pesto_CI, group);
figure(2)
hold on
title('ANOVA order of EC50 PESTO CI - LacIrep3', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_LacIrep3_pesto_CI);
hold on
title('Multcompare order of EC50 PESTO CI - LacIrep3', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep3_table_PESTO_CI', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_Order_LacIrep3_PESTO_CI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC50_Order_LacIrep3_PESTO_CI', '.jpg']);

close all 

