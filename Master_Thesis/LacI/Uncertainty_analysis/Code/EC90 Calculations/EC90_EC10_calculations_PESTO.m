%%Statistics EC90

%SEM
SEM90_1_pesto = std(ec90_1_pesto(:,2))/sqrt(length(ec90_1_pesto(:,2)));
SEM90_2_pesto = std(ec90_2_pesto(:,2))/sqrt(length(ec90_2_pesto(:,2)));
SEM90_3_pesto = std(ec90_3_pesto(:,2))/sqrt(length(ec90_3_pesto(:,2)));

%T-Score
ts90_1_pesto = tinv([0.05 0.95], length(ec90_1_pesto(:,2))-1);
ts90_2_pesto = tinv([0.05 0.95], length(ec90_2_pesto(:,2))-1);
ts90_3_pesto = tinv([0.05 0.95], length(ec90_3_pesto(:,2))-1);

%CI
CI90_1_pesto = mean(ec90_1_pesto(:,2)) + ts90_1_pesto*SEM90_1_pesto;
CI90_2_pesto = mean(ec90_2_pesto(:,2)) + ts90_2_pesto*SEM90_2_pesto;
CI90_3_pesto = mean(ec90_3_pesto(:,2)) + ts90_3_pesto*SEM90_3_pesto;

%STD
std90_1_pesto = std(ec90_1_pesto(:,2));
std90_2_pesto = std(ec90_2_pesto(:,2));
std90_3_pesto = std(ec90_3_pesto(:,2));

%interquantile range
r90_1_pesto = iqr(ec90_1_pesto(:,2));
r90_2_pesto = iqr(ec90_2_pesto(:,2));
r90_3_pesto = iqr(ec90_3_pesto(:,2));

%Mean
mean_ec90_1_pesto = mean(ec90_1_pesto(:,2));
mean_ec90_2_pesto = mean(ec90_2_pesto(:,2));
mean_ec90_3_pesto = mean(ec90_3_pesto(:,2));

%Var
var_ec90_1_pesto = var(ec90_1_pesto(:,2));
var_ec90_2_pesto = var(ec90_2_pesto(:,2));
var_ec90_3_pesto = var(ec90_3_pesto(:,2));

%ANOVA
all_ec90s_pesto = [ec90_1_pesto(:,2) ec90_2_pesto(:,2) ec90_3_pesto(:,2)];
group = ["rep1" "rep2" "rep3"];
[p_90_pesto,tbl_90_pesto,stats_90_pesto] = anova1(all_ec90s_pesto, group);
figure(2)
hold on
title('ANOVA EC90 all Rep - parameter sets from PESTO', 'FontSize',70);
ax = gca;
ax.FontSize = 40;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_all_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_all_PESTO', '.jpg']);
close all

ec90_1_2_pesto = [ec90_1_pesto(:,2) ec90_2_pesto(:,2)];
group = ["rep1" "rep2"];
[p_90_1_2_pesto, tbl_90_1_2_pesto, stats_90_1_2_pesto] = anova1(ec90_1_2_pesto, group);
figure(2)
hold on
title('ANOVA EC90 repression coefficients 1 and 2 - parameter sets from PESTO', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep2_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep2_PESTO', '.jpg']);
close all

ec90_1_3_pesto = [ec90_1_pesto(:,2) ec90_3_pesto(:,2)];
group = ["rep1" "rep3"];
[p_90_1_3_pesto, tbl_90_1_3_pesto, stats_90_1_3_pesto] = anova1(ec90_1_3_pesto, group);
figure(2)
hold on
title('ANOVA EC90 repression coefficients 1 and 3 - parameter sets from PESTO', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep3_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep3_PESTO', '.jpg']);
close all

ec90_2_3_pesto = [ec90_2_pesto(:,2) ec90_3_pesto(:,2)];
group = ["rep2" "rep3"];
[p_90_2_3_pesto, tbl_90_2_3_pesto, stats_90_2_3_pesto] = anova1(ec90_2_3_pesto, group);
figure(2)
hold on
title('ANOVA EC90 repression coefficients 2 and 3 - parameter sets from PESTO', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep2_rep3_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep2_rep3_PESTO', '.jpg']);
close all

%CV
CV90_1_pesto = (std(ec90_1_pesto(:,2))/mean(ec90_1_pesto(:,2))) * 100;
CV90_2_pesto = (std(ec90_2_pesto(:,2))/mean(ec90_2_pesto(:,2))) * 100;
CV90_3_pesto = (std(ec90_3_pesto(:,2))/mean(ec90_3_pesto(:,2))) * 100;

%Histograms
subplot(1,3,1)
hist_ec90_1_pesto = histogram(ec90_1_pesto(:,2));
title('EC90 repression coefficient 1 - parameters from PESTO')
subplot(1,3,2)
hist_ec90_2_pesto = histogram(ec90_2_pesto(:,2));
title('EC90 repression coefficient 2 - parameters from PESTO')
subplot(1,3,3)
hist_ec90_3_pesto = histogram(ec90_3_pesto(:,2));
title('EC90 repression coefficient 3 - parameters from PESTO')

statistics_pesto.sem1 = SEM90_1_pesto;
statistics_pesto.sem2 = SEM90_2_pesto;
statistics_pesto.sem3 = SEM90_3_pesto;

statistics_pesto.tscore1 = ts90_1_pesto;
statistics_pesto.tscore2 = ts90_2_pesto;
statistics_pesto.tscore3 = ts90_3_pesto;

statistics_pesto.CI1 = CI90_1_pesto;
statistics_pesto.CI2 = CI90_2_pesto;
statistics_pesto.CI3 = CI90_3_pesto;

statistics_pesto.std1 = std90_1_pesto;
statistics_pesto.std2 = std90_2_pesto;
statistics_pesto.std3 = std90_3_pesto;

statistics_pesto.interquartile_range1 = r90_1_pesto;
statistics_pesto.interquartile_range2 = r90_2_pesto;
statistics_pesto.interquartile_range3 = r90_3_pesto;

statistics_pesto.mean1 = mean_ec90_1_pesto;
statistics_pesto.mean2 = mean_ec90_2_pesto;
statistics_pesto.mean3 = mean_ec90_3_pesto;

statistics_pesto.var1 = var_ec90_1_pesto;
statistics_pesto.var2 = var_ec90_2_pesto;
statistics_pesto.var3 = var_ec90_3_pesto;

statistics_pesto.ANOVA.all.p = p_90_pesto;
statistics_pesto.ANOVA.all.tbl = tbl_90_pesto;
statistics_pesto.ANOVA.all.stats = stats_90_pesto;

statistics_pesto.ANOVA.rep1_rep3.p = p_90_1_3_pesto;
statistics_pesto.ANOVA.rep1_rep3.tbl = tbl_90_1_3_pesto;
statistics_pesto.ANOVA.rep1_rep3.stats = stats_90_1_3_pesto;

statistics_pesto.ANOVA.rep1_rep2.p = p_90_1_2_pesto;
statistics_pesto.ANOVA.rep1_rep2.tbl = tbl_90_1_2_pesto;
statistics_pesto.ANOVA.rep1_rep2.stats = stats_90_1_2_pesto;

statistics_pesto.ANOVA.rep2_rep3.p = p_90_2_3_pesto;
statistics_pesto.ANOVA.rep2_rep3.tbl = tbl_90_2_3_pesto;
statistics_pesto.ANOVA.rep2_rep3.stats = stats_90_2_3_pesto;

statistics_pesto.CV1 = CV90_1_pesto;
statistics_pesto.CV2 = CV90_2_pesto;
statistics_pesto.CV3 = CV90_3_pesto;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/',datestr(now, 'dd-mmm-yyyy'),'_statistics__EC90_pesto.mat'],'statistics_pesto')

%%Statistics EC10

%SEM
SEM10_1_pesto = std(ec10_1_pesto(:,2))/sqrt(length(ec10_1_pesto(:,2)));
SEM10_2_pesto = std(ec10_2_pesto(:,2))/sqrt(length(ec10_2_pesto(:,2)));
SEM10_3_pesto = std(ec10_3_pesto(:,2))/sqrt(length(ec10_3_pesto(:,2)));

%T-Score
ts10_1_pesto = tinv([0.05 0.95], length(ec10_1_pesto(:,2))-1);
ts10_2_pesto = tinv([0.05 0.95], length(ec10_2_pesto(:,2))-1);
ts10_3_pesto = tinv([0.05 0.95], length(ec10_3_pesto(:,2))-1);

%CI
CI10_1_pesto = mean(ec10_1_pesto(:,2)) + ts10_1_pesto*SEM10_1_pesto;
CI10_2_pesto = mean(ec10_2_pesto(:,2)) + ts10_2_pesto*SEM10_2_pesto;
CI10_3_pesto = mean(ec10_3_pesto(:,2)) + ts10_3_pesto*SEM10_3_pesto;

%STD
std10_1_pesto = std(ec10_1_pesto(:,2));
std10_2_pesto = std(ec10_2_pesto(:,2));
std10_3_pesto = std(ec10_3_pesto(:,2));

%interquantile range
r10_1_pesto = iqr(ec10_1_pesto(:,2));
r10_2_pesto = iqr(ec10_2_pesto(:,2));
r10_3_pesto = iqr(ec10_3_pesto(:,2));

%Mean
mean_ec10_1_pesto = mean(ec10_1_pesto(:,2));
mean_ec10_2_pesto = mean(ec10_2_pesto(:,2));
mean_ec10_3_pesto = mean(ec10_3_pesto(:,2));

%Var
var_ec10_1_pesto = var(ec10_1_pesto(:,2));
var_ec10_2_pesto = var(ec10_2_pesto(:,2));
var_ec10_3_pesto = var(ec10_3_pesto(:,2));

%ANOVA
all_ec10s_pesto = [ec10_1_pesto(:,2) ec10_2_pesto(:,2) ec10_3_pesto(:,2)];
group = ["rep1" "rep2" "rep3"];
[p_10_pesto,tbl_10_pesto,stats_10_pesto] = anova1(all_ec10s_pesto, group);
figure(2)
hold on
title('ANOVA EC10 all Rep - parameter sets from PESTO', 'FontSize',70);
ax = gca;
ax.FontSize = 40;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_all_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_all_PESTO', '.jpg']);
close all

ec10_1_2_pesto = [ec10_1_pesto(:,2) ec10_2_pesto(:,2)];
group = ["rep1" "rep2"];
[p_10_1_2_pesto, tbl_10_1_2_pesto, stats_10_1_2_pesto] = anova1(ec10_1_2_pesto, group);
figure(2)
hold on
title('ANOVA EC10 repression coefficients 1 and 2 - parameter sets from PESTO', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep2_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep2_PESTO', '.jpg']);
close all

ec10_1_3_pesto = [ec10_1_pesto(:,2) ec10_3_pesto(:,2)];
group = ["rep1" "rep3"];
[p_10_1_3_pesto, tbl_10_1_3_pesto, stats_10_1_3_pesto] = anova1(ec10_1_3_pesto, group);
figure(2)
hold on
title('ANOVA EC10 repression coefficients 1 and 3 - parameter sets from PESTO', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep3_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep1_rep3_PESTO', '.jpg']);
close all

ec10_2_3_pesto = [ec10_2_pesto(:,2) ec10_3_pesto(:,2)];
group = ["rep2" "rep3"];
[p_10_2_3_pesto, tbl_10_2_3_pesto, stats_10_2_3_pesto] = anova1(ec10_2_3_pesto, group);
figure(2)
hold on
title('ANOVA EC10 repression coefficients 2 and 3 - parameter sets from PESTO', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep2_rep3_table_PESTO', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC10_rep2_rep3_PESTO', '.jpg']);
close all
%CV
CV10_1_pesto = (std(ec10_1_pesto(:,2))/mean(ec10_1_pesto(:,2))) * 100;
CV10_2_pesto = (std(ec10_2_pesto(:,2))/mean(ec10_2_pesto(:,2))) * 100;
CV10_3_pesto = (std(ec10_3_pesto(:,2))/mean(ec10_3_pesto(:,2))) * 100;

%Histograms
subplot(1,3,1)
hist_ec10_1_pesto = histogram(ec10_1_pesto(:,2));
title('EC10 repression coefficient 1 - parameters from PESTO')
subplot(1,3,2)
hist_ec10_2_pesto = histogram(ec10_2_pesto(:,2));
title('EC10 repression coefficient 2 - parameters from PESTO')
subplot(1,3,3)
hist_ec10_3_pesto = histogram(ec10_3_pesto(:,2));
title('EC10 repression coefficient 3 - parameters from PESTO')

statistics_pesto.sem1 = SEM10_1_pesto;
statistics_pesto.sem2 = SEM10_2_pesto;
statistics_pesto.sem3 = SEM10_3_pesto;

statistics_pesto.tscore1 = ts10_1_pesto;
statistics_pesto.tscore2 = ts10_2_pesto;
statistics_pesto.tscore3 = ts10_3_pesto;

statistics_pesto.CI1 = CI10_1_pesto;
statistics_pesto.CI2 = CI10_2_pesto;
statistics_pesto.CI3 = CI10_3_pesto;

statistics_pesto.std1 = std10_1_pesto;
statistics_pesto.std2 = std10_2_pesto;
statistics_pesto.std3 = std10_3_pesto;

statistics_pesto.interquartile_range1 = r10_1_pesto;
statistics_pesto.interquartile_range2 = r10_2_pesto;
statistics_pesto.interquartile_range3 = r10_3_pesto;

statistics_pesto.mean1 = mean_ec10_1_pesto;
statistics_pesto.mean2 = mean_ec10_2_pesto;
statistics_pesto.mean3 = mean_ec10_3_pesto;

statistics_pesto.var1 = var_ec10_1_pesto;
statistics_pesto.var2 = var_ec10_2_pesto;
statistics_pesto.var3 = var_ec10_3_pesto;

statistics_pesto.ANOVA.all.p = p_10_pesto;
statistics_pesto.ANOVA.all.tbl = tbl_10_pesto;
statistics_pesto.ANOVA.all.stats = stats_10_pesto;

statistics_pesto.ANOVA.rep1_rep3.p = p_10_1_3_pesto;
statistics_pesto.ANOVA.rep1_rep3.tbl = tbl_10_1_3_pesto;
statistics_pesto.ANOVA.rep1_rep3.stats = stats_10_1_3_pesto;

statistics_pesto.ANOVA.rep1_rep2.p = p_10_1_2_pesto;
statistics_pesto.ANOVA.rep1_rep2.tbl = tbl_10_1_2_pesto;
statistics_pesto.ANOVA.rep1_rep2.stats = stats_10_1_2_pesto;

statistics_pesto.ANOVA.rep2_rep3.p = p_10_2_3_pesto;
statistics_pesto.ANOVA.rep2_rep3.tbl = tbl_10_2_3_pesto;
statistics_pesto.ANOVA.rep2_rep3.stats = stats_10_2_3_pesto;

statistics_pesto.CV1 = CV10_1_pesto;
statistics_pesto.CV2 = CV10_2_pesto;
statistics_pesto.CV3 = CV10_3_pesto;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/',datestr(now, 'dd-mmm-yyyy'),'_statistics__EC10_PESTO.mat'],'statistics_pesto')
