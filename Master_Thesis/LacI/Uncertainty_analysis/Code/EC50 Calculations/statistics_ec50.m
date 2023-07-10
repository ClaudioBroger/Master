%standard error
SEM50_1 = std(ec50_1(:,2))/sqrt(length(ec50_1(:,2)));
SEM50_2 = std(ec50_2(:,2))/sqrt(length(ec50_2(:,2)));
SEM50_3 = std(ec50_3(:,2))/sqrt(length(ec50_3(:,2)));

SEM90_1 = std(ec90_1(:,2))/sqrt(length(ec90_1(:,2)));
SEM90_2 = std(ec90_2(:,2))/sqrt(length(ec90_2(:,2)));
SEM90_3 = std(ec90_3(:,2))/sqrt(length(ec90_3(:,2)));

%T-Score
ts50_1 = tinv([0.05 0.95], length(ec50_1(:,2))-1);
ts50_2 = tinv([0.05 0.95], length(ec50_2(:,2))-1);
ts50_3 = tinv([0.05 0.95], length(ec50_3(:,2))-1);

ts90_1 = tinv([0.05 0.95], length(ec90_1(:,2))-1);
ts90_2 = tinv([0.05 0.95], length(ec90_2(:,2))-1);
ts90_3 = tinv([0.05 0.95], length(ec90_3(:,2))-1);

%Confidence Interval
CI50_1 = mean(ec50_1(:,2)) + ts50_1*SEM50_1;
CI50_2 = mean(ec50_2(:,2)) + ts50_2*SEM50_2;
CI50_3 = mean(ec50_3(:,2)) + ts50_3*SEM50_3;

CI90_1 = mean(ec90_1(:,2)) + ts90_1*SEM90_1;
CI90_2 = mean(ec90_2(:,2)) + ts90_2*SEM90_2;
CI90_3 = mean(ec90_3(:,2)) + ts90_3*SEM90_3;

%Standard deviation
std50_1 = std(ec50_1(:,2));
std50_2 = std(ec50_2(:,2));
std50_3 = std(ec50_3(:,2));

std90_1 = std(ec90_1(:,2));
std90_2 = std(ec90_2(:,2));
std90_3 = std(ec90_3(:,2));

%Interquantile range
r50_1 = iqr(ec50_1(:,2));
r50_2 = iqr(ec50_2(:,2));
r50_3 = iqr(ec50_3(:,2));

r90_1 = iqr(ec90_1(:,2));
r90_2 = iqr(ec90_2(:,2));
r90_3 = iqr(ec90_3(:,2));

%Mean
mean_ec50_1 = mean(ec50_1(:,2));
mean_ec50_2 = mean(ec50_2(:,2));
mean_ec50_3 = mean(ec50_3(:,2));

mean_ec90_1 = mean(ec90_1(:,2));
mean_ec90_2 = mean(ec90_2(:,2));
mean_ec90_3 = mean(ec90_3(:,2));

%Variance
var_ec50_1 = var(ec50_1(:,2));
var_ec50_2 = var(ec50_2(:,2));
var_ec50_3 = var(ec50_3(:,2));

var_ec90_1 = var(ec90_1(:,2));
var_ec90_2 = var(ec90_2(:,2));
var_ec90_3 = var(ec90_3(:,2));


%ANOVA
all_ec50s = [ec50_1(:,2) ec50_2(:,2) ec50_3(:,2)];
group = ["rep1" "rep2" "rep3"];
[p_50,tbl_50,stats_50] = anova1(all_ec50s, group);
figure(2)
hold on
title('ANOVA EC50 all Rep - parameter sets from hyperspace', 'FontSize',70);
ax = gca;
ax.FontSize = 40;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_alltable', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_all', '.jpg']);
close all

ec50_1_2 = [ec50_1(:,2) ec50_2(:,2)];
group = ["rep1" "rep2"];
[p_50_1_2, tbl_50_1_2, stats_50_1_2] = anova1(ec50_1_2, group);
figure(2)
hold on
title('ANOVA EC50 repression coefficients 1 and 2 - parameter sets from hyperspace', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_rep1_rep2_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_rep1_rep2', '.jpg']);
close all

ec50_1_3 = [ec50_1(:,2) ec50_3(:,2)];
group = ["rep1" "rep3"];
[p_50_1_3, tbl_50_1_3, stats_50_1_3] = anova1(ec50_1_3, group);
figure(2)
hold on
title('ANOVA EC50 repression coefficients 1 and 3 - parameter sets from hyperspace', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_rep1_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_rep1_rep3', '.jpg']);
close all

ec50_2_3 = [ec50_2(:,2) ec50_3(:,2)];
group = ["rep2" "rep3"];
[p_50_2_3, tbl_50_2_3, stats_50_2_3] = anova1(ec50_2_3, group);
figure(2)
hold on
title('ANOVA EC50 repression coefficients 2 and 3 - parameter sets from hyperspace', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_rep2_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC50_rep2_rep3', '.jpg']);
close all


%Coefficient of Variation
CV50_1 = (std(ec50_1(:,2))/mean(ec50_1(:,2))) * 100;
CV50_2 = (std(ec50_2(:,2))/mean(ec50_2(:,2))) * 100;
CV50_3 = (std(ec50_3(:,2))/mean(ec50_3(:,2))) * 100;

CV90_1 = (std(ec90_1(:,2))/mean(ec90_1(:,2))) * 100;
CV90_2 = (std(ec90_2(:,2))/mean(ec90_2(:,2))) * 100;
CV90_3 = (std(ec90_3(:,2))/mean(ec90_3(:,2))) * 100;

%Histograms
subplot(1,3,1)
hist_ec50_1 = histogram(ec50_1(:,2));
title('EC50 repression coefficient 1 - randomly drawn parameters')
subplot(1,3,2)
hist_ec50_2 = histogram(ec50_2(:,2));
title('EC50 repression coefficient 2 - randomly drawn parameters')
subplot(1,3,3)
hist_ec50_3 = histogram(ec50_3(:,2));
title('EC50 repression coefficient 3 - randomly drawn parameters')
subplot(2,3,4)
hist_ec50_1_around_paraopt = histogram(ec50_1_around_paraopt(:,2));
title('EC50 repression coefficient 1 - parameters around optimal parameterset')
subplot(2,3,5)
hist_ec50_2_around_paraopt = histogram(ec50_2_around_paraopt(:,2));
title('EC50 repression coefficient 2 - parameters around optimal parameterset')
subplot(2,3,6)
hist_ec50_3_around_paraopt = histogram(ec50_3_around_paraopt(:,2));
title('EC50 repression coefficient 3 - parameters around optimal parameterset')


subplot(1,3,1)
hist_ec90_1 = histogram(ec90_1(:,2));
title('EC90 repression coefficient 1 - randomly drawn parameters')
subplot(1,3,2)
hist_ec90_2 = histogram(ec90_2(:,2));
title('EC90 repression coefficient 2 - randomly drawn parameters')
subplot(1,3,3)
hist_ec90_3 = histogram(ec90_3(:,2));
title('EC90 repression coefficient 3 - randomly drawn parameters')


%Check if the order of dose response curves for the same parameters but
%with different repression coefficient is the same all the time.

for k = 1:length(ec50_1(:,2))
    if (ec50_1(k,2) < ec50_2(k,2)) && (ec50_2(k,2)< ec50_3(k,2))
        order(k) = 1;
    else
        order(k) = 0;
    end
end

%save in one structure

statistics_hyperspace.sem1 = SEM50_1;
statistics_hyperspace.sem2 = SEM50_2;
statistics_hyperspace.sem3 = SEM50_3;

statistics_hyperspace.tscore1 = ts50_1;
statistics_hyperspace.tscore2 = ts50_2;
statistics_hyperspace.tscore3 = ts50_3;

statistics_hyperspace.CI1 = CI50_1;
statistics_hyperspace.CI2 = CI50_2;
statistics_hyperspace.CI3 = CI50_3;

statistics_hyperspace.std1 = std50_1;
statistics_hyperspace.std2 = std50_2;
statistics_hyperspace.std3 = std50_3;

statistics_hyperspace.interquartile_range1 = r50_1;
statistics_hyperspace.interquartile_range2 = r50_2;
statistics_hyperspace.interquartile_range3 = r50_3;

statistics_hyperspace.mean1 = mean_ec50_1;
statistics_hyperspace.mean2 = mean_ec50_2;
statistics_hyperspace.mean3 = mean_ec50_3;

statistics_hyperspace.var1 = var_ec50_1;
statistics_hyperspace.var2 = var_ec50_2;
statistics_hyperspace.var3 = var_ec50_3;

statistics_hyperspace.ANOVA.all.p = p_50;
statistics_hyperspace.ANOVA.all.tbl = tbl_50;
statistics_hyperspace.ANOVA.all.stats = stats_50;

statistics_hyperspace.ANOVA.rep1_rep3.p = p_50_1_3;
statistics_hyperspace.ANOVA.rep1_rep3.tbl = tbl_50_1_3;
statistics_hyperspace.ANOVA.rep1_rep3.stats = stats_50_1_3;

statistics_hyperspace.ANOVA.rep1_rep2.p = p_50_1_2;
statistics_hyperspace.ANOVA.rep1_rep2.tbl = tbl_50_1_2;
statistics_hyperspace.ANOVA.rep1_rep2.stats = stats_50_1_2;

statistics_hyperspace.ANOVA.rep2_rep3.p = p_50_2_3;
statistics_hyperspace.ANOVA.rep2_rep3.tbl = tbl_50_2_3;
statistics_hyperspace.ANOVA.rep2_rep3.stats = stats_50_2_3;

statistics_hyperspace.CV1 = CV50_1;
statistics_hyperspace.CV2 = CV50_2;
statistics_hyperspace.CV3 = CV50_3;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from Hyperspace/',datestr(now, 'dd-mmm-yyyy'),'_statistics_hyperspace.mat'],'statistics_hyperspace')




