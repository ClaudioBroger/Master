%% Calculate order of EC50 and differences in parameters 

%%Hyperspace

load('28-Apr-2023ec10_1_TetRTup1.mat');
load('28-Apr-2023ec10_2_TetRTup1.mat');
load('29-Apr-2023ec10_3_TetRTup1.mat');

ec10_1over2 = find(ec10_1(:,2) >= ec10_2(:,2));
ec10_2over3 = find(ec10_2(:,2) >= ec10_3(:,2));
ec10_1over3 = find(ec10_1(:,2) >= ec10_3(:,2));


%9 EC10 found rep 1 over rep 2


%15 EC50 from first repression coefficient found to be over EC50 from
%second repression coefficient

load('27-Apr-2023_rand_parameter_LacI_TetR_EC50.mat');

rand_parameter_wrong_order = rand_parameter(ec10_1over2,:);

rand_parameter_right_order = rand_parameter;
rand_parameter_right_order(ec10_1over2,:) = [];

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

kLacI = [table2array(rand_parameter_right_order(:,"kLacI")) table2array(rand_parameter_wrong_order(:,"kLacI"))];
group = ["right" "wrong"];
[p_kLacI, tbl_kLacI, stats_kLacI] = anova1(kLacI, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - kLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_kLacI);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - kLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_kLacI_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_kLacI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_kLacI', '.jpg']);

close all 


kCit = [table2array(rand_parameter_right_order(:,"kCit")) table2array(rand_parameter_wrong_order(:,"kCit"))];
group = ["right" "wrong"];
[p_kCit, tbl_kCit, stats_kCit] = anova1(kCit, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - kCit', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_kCit);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - kCit', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_kCit_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_kCit', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_kCit', '.jpg']);
close all  

dLacI = [table2array(rand_parameter_right_order(:,"dLacI")) table2array(rand_parameter_wrong_order(:,"dLacI"))];
group = ["right" "wrong"];
[p_dLacI, tbl_dLacI, stats_dLacI] = anova1(dLacI, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - dLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_dLacI);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - dLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_dLacI_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_dLacI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_dLacI', '.jpg']);

close all  

LacIrep = [table2array(rand_parameter_right_order(:,"LacIrep")) table2array(rand_parameter_wrong_order(:,"LacIrep"))];
group = ["right" "wrong"];
[p_LacIrep, tbl_LacIrep, stats_LacIrep] = anova1(LacIrep, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - LacIrep', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_LacIrep);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - LacIrep', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_LacIrep_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_LacIrep', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_LacIrep', '.jpg']);

close all  

nLacI = [table2array(rand_parameter_right_order(:,"nLacI")) table2array(rand_parameter_wrong_order(:,"nLacI"))];
group = ["right" "wrong"];
[p_nLacI, tbl_nLacI, stats_nLacI] = anova1(nLacI, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - nLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_nLacI);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - nLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_nLacI_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_nLacI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_nLacI', '.jpg']);

close all  

CitL = [table2array(rand_parameter_right_order(:,"CitL")) table2array(rand_parameter_wrong_order(:,"CitL"))];
group = ["right" "wrong"];
[p_CitL, tbl_CitL, stats_CitL] = anova1(CitL, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - CitL', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_CitL);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - CitL', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_CitL_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_CitL', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_CitL', '.jpg']);

close all  

KdLacI = [table2array(rand_parameter_right_order(:,"KdLacI")) table2array(rand_parameter_wrong_order(:,"KdLacI"))];
group = ["right" "wrong"];
[p_KdLacI, tbl_KdLacI, stats_KdLacI] = anova1(KdLacI, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - KdLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_KdLacI);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - KdLacI', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_KdLacI_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_KdLacI', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_KdLacI', '.jpg']);

close all 

LacIrep2 = [table2array(rand_parameter_right_order(:,"LacIrep2")) table2array(rand_parameter_wrong_order(:,"LacIrep2"))];
group = ["right" "wrong"];
[p_LacIrep2, tbl_LacIrep2, stats_LacIrep2] = anova1(LacIrep2, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - LacIrep2', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_LacIrep2);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - LacIrep2', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_LacIrep2_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_LacIrep2', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_LacIrep2', '.jpg']);

close all 

dTetRTup1 = [table2array(rand_parameter_right_order(:,"LacIrep3")) table2array(rand_parameter_wrong_order(:,"LacIrep3"))];
group = ["right" "wrong"];
[p_LacIrep3, tbl_LacIrep3, stats_LacIrep3] = anova1(dTetRTup1, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - LacIrep3', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_LacIrep3);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - LacIrep3', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_LacIrep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_LacIrep3', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_LacIrep3', '.jpg']);

close all 


kTetRTup1 = [table2array(rand_parameter_right_order(:,"kTetRTup1")) table2array(rand_parameter_wrong_order(:,"kTetRTup1"))];
group = ["right" "wrong"];
[p_kTetRTup1, tbl_kTetRTup1, stats_kTetRTup1] = anova1(kTetRTup1, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - kTetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_kTetRTup1);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - kTetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_kTetRTup1_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_kTetRTup1', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_kTetRTup1', '.jpg']);

close all 


dTetRTup1 = [table2array(rand_parameter_right_order(:,"dTetRTup1")) table2array(rand_parameter_wrong_order(:,"dTetRTup1"))];
group = ["right" "wrong"];
[p_dTetRTup1, tbl_dTetRTup1, stats_dTetRTup1] = anova1(dTetRTup1, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - dTetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_dTetRTup1);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - dTetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_dTetRTup1_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_dTetRTup1', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_dTetRTup1', '.jpg']);

close all 


degtag = [table2array(rand_parameter_right_order(:,"degtag")) table2array(rand_parameter_wrong_order(:,"degtag"))];
group = ["right" "wrong"];
[p_degtag, tbl_degtag, stats_degtag] = anova1(degtag, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - degtag', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_degtag);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - degtag', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_degtag_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_degtag', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_degtag', '.jpg']);

close all 


TetRTup1L = [table2array(rand_parameter_right_order(:,"TetRTup1L")) table2array(rand_parameter_wrong_order(:,"TetRTup1L"))];
group = ["right" "wrong"];
[p_TetRTup1L, tbl_TetRTup1L, stats_TetRTup1L] = anova1(TetRTup1L, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - TetRTup1L', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_TetRTup1L);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - TetRTup1L', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_TetRTup1L_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_TetRTup1L', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_TetRTup1L', '.jpg']);

close all 


nTetRTup1 = [table2array(rand_parameter_right_order(:,"nTetRTup1")) table2array(rand_parameter_wrong_order(:,"nTetRTup1"))];
group = ["right" "wrong"];
[p_nTetRTup1, tbl_nTetRTup1, stats_nTetRTup1] = anova1(nTetRTup1, group);
figure(2)
hold on
title('ANOVA order of EC10 LacI & TetRTup1 - nTetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
figure(3)
c = multcompare(stats_nTetRTup1);
hold on
title('Multcompare order of EC10 LacI & TetRTup1 - nTetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
set(figure(3), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_nTetRTup1_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_Order_nTetRTup1', '.jpg']);
saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'Multcompare_EC10_Order_nTetRTup1', '.jpg']);

close all 



