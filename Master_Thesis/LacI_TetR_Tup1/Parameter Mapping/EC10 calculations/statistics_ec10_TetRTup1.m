%standard error
SEM10_1 = std(ec10_1(:,2))/sqrt(length(ec10_1(:,2)));
SEM10_2 = std(ec10_2(:,2))/sqrt(length(ec10_2(:,2)));
SEM10_3 = std(ec10_3(:,2))/sqrt(length(ec10_3(:,2)));


%T-Score
ts10_1 = tinv([0.05 0.95], length(ec10_1(:,2))-1);
ts10_2 = tinv([0.05 0.95], length(ec10_2(:,2))-1);
ts10_3 = tinv([0.05 0.95], length(ec10_3(:,2))-1);


%Confidence Interval
CI10_1 = mean(ec10_1(:,2)) + ts10_1*SEM10_1;
CI10_2 = mean(ec10_2(:,2)) + ts10_2*SEM10_2;
CI10_3 = mean(ec10_3(:,2)) + ts10_3*SEM10_3;

%Standard deviation
std10_1 = std(ec10_1(:,2));
std10_2 = std(ec10_2(:,2));
std10_3 = std(ec10_3(:,2));


%Interquantile range
r10_1 = iqr(ec10_1(:,2));
r10_2 = iqr(ec10_2(:,2));
r10_3 = iqr(ec10_3(:,2));


%Mean
mean_ec10_1 = mean(ec10_1(:,2));
mean_ec10_2 = mean(ec10_2(:,2));
mean_ec10_3 = mean(ec10_3(:,2));

%Variance
var_ec10_1 = var(ec10_1(:,2));
var_ec10_2 = var(ec10_2(:,2));
var_ec10_3 = var(ec10_3(:,2));


%ANOVA
all_ec10s = [ec10_1(:,2) ec10_2(:,2) ec10_3(:,2)];
group = ["rep1" "rep2" "rep3"];
[p_10,tbl_10,stats_10] = anova1(all_ec10s, group);
figure(2)
hold on
title('ANOVA all repression coefficients - LacI & TetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_alltable', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_all', '.jpg']);
close all

ec10_1_2 = [ec10_1(:,2) ec10_2(:,2)];
group = ["rep1" "rep2"];
[p_10_1_2, tbl_10_1_2, stats_10_1_2] = anova1(ec10_1_2, group);
figure(2)
hold on
title('ANOVA repression coefficients 1 and 2 - LacI & TetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep2_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep2', '.jpg']);
close all

ec10_1_3 = [ec10_1(:,2) ec10_3(:,2)];
group = ["rep1" "rep3"];
[p_10_1_3, tbl_10_1_3, stats_10_1_3] = anova1(ec10_1_3, group);
figure(2)
hold on
title('ANOVA repression coefficients 1 and 3 - parameter sets from hyperspace', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep3', '.jpg']);
close all

ec10_2_3 = [ec10_2(:,2) ec10_3(:,2)];
group = ["rep2" "rep3"];
[p_10_2_3, tbl_10_2_3, stats_10_2_3] = anova1(ec10_2_3, group);
figure(2)
hold on
title('ANOVA repression coefficients 2 and 3 - LacI & TetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep2_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep2_rep3', '.jpg']);
close all


%Coefficient of Variation
CV10_1 = (std(ec10_1(:,2))/mean(ec10_1(:,2))) * 100;
CV10_2 = (std(ec10_2(:,2))/mean(ec10_2(:,2))) * 100;
CV10_3 = (std(ec10_3(:,2))/mean(ec10_3(:,2))) * 100;


%Histograms
subplot(1,3,1)
hist_ec10_1 = histogram(ec10_1(:,2));
title('EC10 repression coefficient 1 - randomly drawn parameters')
subplot(1,3,2)
hist_ec10_2 = histogram(ec10_2(:,2));
title('EC10 repression coefficient 2 - randomly drawn parameters')
subplot(1,3,3)
hist_ec10_3 = histogram(ec10_3(:,2));
title('EC10 repression coefficient 3 - randomly drawn parameters')







%save in one structure

statistics_TetRTup1_EC10.sem1 = SEM10_1;
statistics_TetRTup1_EC10.sem2 = SEM10_2;
statistics_TetRTup1_EC10.sem3 = SEM10_3;

statistics_TetRTup1_EC10.tscore1 = ts10_1;
statistics_TetRTup1_EC10.tscore2 = ts10_2;
statistics_TetRTup1_EC10.tscore3 = ts10_3;

statistics_TetRTup1_EC10.CI1 = CI10_1;
statistics_TetRTup1_EC10.CI2 = CI10_2;
statistics_TetRTup1_EC10.CI3 = CI10_3;

statistics_TetRTup1_EC10.std1 = std10_1;
statistics_TetRTup1_EC10.std2 = std10_2;
statistics_TetRTup1_EC10.std3 = std10_3;

statistics_TetRTup1_EC10.interquartile_range1 = r10_1;
statistics_TetRTup1_EC10.interquartile_range2 = r10_2;
statistics_TetRTup1_EC10.interquartile_range3 = r10_3;

statistics_TetRTup1_EC10.mean1 = mean_ec10_1;
statistics_TetRTup1_EC10.mean2 = mean_ec10_2;
statistics_TetRTup1_EC10.mean3 = mean_ec10_3;

statistics_TetRTup1_EC10.var1 = var_ec10_1;
statistics_TetRTup1_EC10.var2 = var_ec10_2;
statistics_TetRTup1_EC10.var3 = var_ec10_3;

statistics_TetRTup1_EC10.ANOVA.all.p = p_10;
statistics_TetRTup1_EC10.ANOVA.all.tbl = tbl_10;
statistics_TetRTup1_EC10.ANOVA.all.stats = stats_10;

statistics_TetRTup1_EC10.ANOVA.rep1_rep3.p = p_10_1_3;
statistics_TetRTup1_EC10.ANOVA.rep1_rep3.tbl = tbl_10_1_3;
statistics_TetRTup1_EC10.ANOVA.rep1_rep3.stats = stats_10_1_3;

statistics_TetRTup1_EC10.ANOVA.rep1_rep2.p = p_10_1_2;
statistics_TetRTup1_EC10.ANOVA.rep1_rep2.tbl = tbl_10_1_2;
statistics_TetRTup1_EC10.ANOVA.rep1_rep2.stats = stats_10_1_2;

statistics_TetRTup1_EC10.ANOVA.rep2_rep3.p = p_10_2_3;
statistics_TetRTup1_EC10.ANOVA.rep2_rep3.tbl = tbl_10_2_3;
statistics_TetRTup1_EC10.ANOVA.rep2_rep3.stats = stats_10_2_3;

statistics_TetRTup1_EC10.CV1 = CV10_1;
statistics_TetRTup1_EC10.CV2 = CV10_2;
statistics_TetRTup1_EC10.CV3 = CV10_3;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/',datestr(now, 'dd-mmm-yyyy'),'_statistics_TetRTup1_EC10.mat'],'statistics_TetRTup1_EC10')




