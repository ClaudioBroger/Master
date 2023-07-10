%%Statistics EC90

%SEM
SEM90_1 = std(ec90_1(:,2))/sqrt(length(ec90_1(:,2)));
SEM90_2 = std(ec90_2(:,2))/sqrt(length(ec90_2(:,2)));
SEM90_3 = std(ec90_3(:,2))/sqrt(length(ec90_3(:,2)));

%T-Score
ts90_1 = tinv([0.05 0.95], length(ec90_1(:,2))-1);
ts90_2 = tinv([0.05 0.95], length(ec90_2(:,2))-1);
ts90_3 = tinv([0.05 0.95], length(ec90_3(:,2))-1);

%CI
CI90_1 = mean(ec90_1(:,2)) + ts90_1*SEM90_1;
CI90_2 = mean(ec90_2(:,2)) + ts90_2*SEM90_2;
CI90_3 = mean(ec90_3(:,2)) + ts90_3*SEM90_3;

%STD
std90_1 = std(ec90_1(:,2));
std90_2 = std(ec90_2(:,2));
std90_3 = std(ec90_3(:,2));

%interquantile range
r90_1 = iqr(ec90_1(:,2));
r90_2 = iqr(ec90_2(:,2));
r90_3 = iqr(ec90_3(:,2));

%Mean
mean_ec90_1 = mean(ec90_1(:,2));
mean_ec90_2 = mean(ec90_2(:,2));
mean_ec90_3 = mean(ec90_3(:,2));

%Var
var_ec90_1 = var(ec90_1(:,2));
var_ec90_2 = var(ec90_2(:,2));
var_ec90_3 = var(ec90_3(:,2));

%ANOVA
all_ec90s = [ec90_1(:,2) ec90_2(:,2) ec90_3(:,2)];
group = ["rep1" "rep2" "rep3"];
[p_90,tbl_90,stats_90] = anova1(all_ec90s, group);
figure(2)
hold on
title('ANOVA EC90 all Rep - parameter sets from hyperspace', 'FontSize',70);
ax = gca;
ax.FontSize = 40;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_alltable', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_all', '.jpg']);
close all

ec90_1_2 = [ec90_1(:,2) ec90_2(:,2)];
group = ["rep1" "rep2"];
[p_90_1_2, tbl_90_1_2, stats_90_1_2] = anova1(ec90_1_2, group);
figure(2)
hold on
title('ANOVA EC90 repression coefficients 1 and 2 - parameter sets from hyperspace', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep2_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep2', '.jpg']);
close all

ec90_1_3 = [ec90_1(:,2) ec90_3(:,2)];
group = ["rep1" "rep3"];
[p_90_1_3, tbl_90_1_3, stats_90_1_3] = anova1(ec90_1_3, group);
figure(2)
hold on
title('ANOVA EC90 repression coefficients 1 and 3 - parameter sets from hyperspace', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep3', '.jpg']);
close all

ec90_2_3 = [ec90_2(:,2) ec90_3(:,2)];
group = ["rep2" "rep3"];
[p_90_2_3, tbl_90_2_3, stats_90_2_3] = anova1(ec90_2_3, group);
figure(2)
hold on
title('ANOVA EC90 repression coefficients 2 and 3 - parameter sets from hyperspace', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep2_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep2_rep3', '.jpg']);
close all

%CV
CV90_1 = (std(ec90_1(:,2))/mean(ec90_1(:,2))) * 100;
CV90_2 = (std(ec90_2(:,2))/mean(ec90_2(:,2))) * 100;
CV90_3 = (std(ec90_3(:,2))/mean(ec90_3(:,2))) * 100;

%Histograms
subplot(1,3,1)
hist_ec90_1 = histogram(ec90_1(:,2));
title('EC90 repression coefficient 1 - randomly drawn parameters')
subplot(1,3,2)
hist_ec90_2 = histogram(ec90_2(:,2));
title('EC90 repression coefficient 2 - randomly drawn parameters')
subplot(1,3,3)
hist_ec90_3 = histogram(ec90_3(:,2));
title('EC90 repression coefficient 3 - randomly drawn parameters')

statistics_hyperspace.sem1 = SEM90_1;
statistics_hyperspace.sem2 = SEM90_2;
statistics_hyperspace.sem3 = SEM90_3;

statistics_hyperspace.tscore1 = ts90_1;
statistics_hyperspace.tscore2 = ts90_2;
statistics_hyperspace.tscore3 = ts90_3;

statistics_hyperspace.CI1 = CI90_1;
statistics_hyperspace.CI2 = CI90_2;
statistics_hyperspace.CI3 = CI90_3;

statistics_hyperspace.std1 = std90_1;
statistics_hyperspace.std2 = std90_2;
statistics_hyperspace.std3 = std90_3;

statistics_hyperspace.interquartile_range1 = r90_1;
statistics_hyperspace.interquartile_range2 = r90_2;
statistics_hyperspace.interquartile_range3 = r90_3;

statistics_hyperspace.mean1 = mean_ec90_1;
statistics_hyperspace.mean2 = mean_ec90_2;
statistics_hyperspace.mean3 = mean_ec90_3;

statistics_hyperspace.var1 = var_ec90_1;
statistics_hyperspace.var2 = var_ec90_2;
statistics_hyperspace.var3 = var_ec90_3;

statistics_hyperspace.ANOVA.all.p = p_90;
statistics_hyperspace.ANOVA.all.tbl = tbl_90;
statistics_hyperspace.ANOVA.all.stats = stats_90;

statistics_hyperspace.ANOVA.rep1_rep3.p = p_90_1_3;
statistics_hyperspace.ANOVA.rep1_rep3.tbl = tbl_90_1_3;
statistics_hyperspace.ANOVA.rep1_rep3.stats = stats_90_1_3;

statistics_hyperspace.ANOVA.rep1_rep2.p = p_90_1_2;
statistics_hyperspace.ANOVA.rep1_rep2.tbl = tbl_90_1_2;
statistics_hyperspace.ANOVA.rep1_rep2.stats = stats_90_1_2;

statistics_hyperspace.ANOVA.rep2_rep3.p = p_90_2_3;
statistics_hyperspace.ANOVA.rep2_rep3.tbl = tbl_90_2_3;
statistics_hyperspace.ANOVA.rep2_rep3.stats = stats_90_2_3;

statistics_hyperspace.CV1 = CV90_1;
statistics_hyperspace.CV2 = CV90_2;
statistics_hyperspace.CV3 = CV90_3;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/',datestr(now, 'dd-mmm-yyyy'),'_statistics__EC90_hyperspace.mat'],'statistics_hyperspace')

%%Statistics EC10

%SEM
SEM10_1 = std(ec10_1(:,2))/sqrt(length(ec10_1(:,2)));
SEM10_2 = std(ec10_2(:,2))/sqrt(length(ec10_2(:,2)));
SEM10_3 = std(ec10_3(:,2))/sqrt(length(ec10_3(:,2)));

%T-Score
ts10_1 = tinv([0.05 0.95], length(ec10_1(:,2))-1);
ts10_2 = tinv([0.05 0.95], length(ec10_2(:,2))-1);
ts10_3 = tinv([0.05 0.95], length(ec10_3(:,2))-1);

%CI
CI10_1 = mean(ec10_1(:,2)) + ts10_1*SEM10_1;
CI10_2 = mean(ec10_2(:,2)) + ts10_2*SEM10_2;
CI10_3 = mean(ec10_3(:,2)) + ts10_3*SEM10_3;

%STD
std10_1 = std(ec10_1(:,2));
std10_2 = std(ec10_2(:,2));
std10_3 = std(ec10_3(:,2));

%interquantile range
r10_1 = iqr(ec10_1(:,2));
r10_2 = iqr(ec10_2(:,2));
r10_3 = iqr(ec10_3(:,2));

%Mean
mean_ec10_1 = mean(ec10_1(:,2));
mean_ec10_2 = mean(ec10_2(:,2));
mean_ec10_3 = mean(ec10_3(:,2));

%Var
var_ec10_1 = var(ec10_1(:,2));
var_ec10_2 = var(ec10_2(:,2));
var_ec10_3 = var(ec10_3(:,2));

%ANOVA
all_ec10s = [ec10_1(:,2) ec10_2(:,2) ec10_3(:,2)];
group = ["rep1" "rep2" "rep3"];
[p_10,tbl_10,stats_10] = anova1(all_ec10s, group);
figure(2)
hold on
title('ANOVA EC10 all Rep - parameter sets from hyperspace', 'FontSize',70);
ax = gca;
ax.FontSize = 40;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_alltable', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_all', '.jpg']);
close all

ec10_1_2 = [ec10_1(:,2) ec10_2(:,2)];
group = ["rep1" "rep2"];
[p_10_1_2, tbl_10_1_2, stats_10_1_2] = anova1(ec10_1_2, group);
figure(2)
hold on
title('ANOVA EC10 repression coefficients 1 and 2 - parameter sets from hyperspace', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep2_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep2', '.jpg']);
close all

ec10_1_3 = [ec10_1(:,2) ec10_3(:,2)];
group = ["rep1" "rep3"];
[p_10_1_3, tbl_10_1_3, stats_10_1_3] = anova1(ec10_1_3, group);
figure(2)
hold on
title('ANOVA EC10 repression coefficients 1 and 3 - parameter sets from hyperspace', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep3', '.jpg']);
close all

ec10_2_3 = [ec10_2(:,2) ec10_3(:,2)];
group = ["rep2" "rep3"];
[p_10_2_3, tbl_10_2_3, stats_10_2_3] = anova1(ec10_2_3, group);
figure(2)
hold on
title('ANOVA EC10 repression coefficients 2 and 3 - parameter sets from hyperspace', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep2_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep2_rep3', '.jpg']);
close all

%CV
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

statistics_hyperspace.sem1 = SEM10_1;
statistics_hyperspace.sem2 = SEM10_2;
statistics_hyperspace.sem3 = SEM10_3;

statistics_hyperspace.tscore1 = ts10_1;
statistics_hyperspace.tscore2 = ts10_2;
statistics_hyperspace.tscore3 = ts10_3;

statistics_hyperspace.CI1 = CI10_1;
statistics_hyperspace.CI2 = CI10_2;
statistics_hyperspace.CI3 = CI10_3;

statistics_hyperspace.std1 = std10_1;
statistics_hyperspace.std2 = std10_2;
statistics_hyperspace.std3 = std10_3;

statistics_hyperspace.interquartile_range1 = r10_1;
statistics_hyperspace.interquartile_range2 = r10_2;
statistics_hyperspace.interquartile_range3 = r10_3;

statistics_hyperspace.mean1 = mean_ec10_1;
statistics_hyperspace.mean2 = mean_ec10_2;
statistics_hyperspace.mean3 = mean_ec10_3;

statistics_hyperspace.var1 = var_ec10_1;
statistics_hyperspace.var2 = var_ec10_2;
statistics_hyperspace.var3 = var_ec10_3;

statistics_hyperspace.ANOVA.all.p = p_10;
statistics_hyperspace.ANOVA.all.tbl = tbl_10;
statistics_hyperspace.ANOVA.all.stats = stats_10;

statistics_hyperspace.ANOVA.rep1_rep3.p = p_10_1_3;
statistics_hyperspace.ANOVA.rep1_rep3.tbl = tbl_10_1_3;
statistics_hyperspace.ANOVA.rep1_rep3.stats = stats_10_1_3;

statistics_hyperspace.ANOVA.rep1_rep2.p = p_10_1_2;
statistics_hyperspace.ANOVA.rep1_rep2.tbl = tbl_10_1_2;
statistics_hyperspace.ANOVA.rep1_rep2.stats = stats_10_1_2;

statistics_hyperspace.ANOVA.rep2_rep3.p = p_10_2_3;
statistics_hyperspace.ANOVA.rep2_rep3.tbl = tbl_10_2_3;
statistics_hyperspace.ANOVA.rep2_rep3.stats = stats_10_2_3;

statistics_hyperspace.CV1 = CV10_1;
statistics_hyperspace.CV2 = CV10_2;
statistics_hyperspace.CV3 = CV10_3;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/',datestr(now, 'dd-mmm-yyyy'),'_statistics__EC10_hyperspace.mat'],'statistics_hyperspace')
