%standard error
SEM90_1 = std(ec90_1(:,2))/sqrt(length(ec90_1(:,2)));
SEM90_2 = std(ec90_2(:,2))/sqrt(length(ec90_2(:,2)));
SEM90_3 = std(ec90_3(:,2))/sqrt(length(ec90_3(:,2)));



%T-Score
ts90_1 = tinv([0.05 0.95], length(ec90_1(:,2))-1);
ts90_2 = tinv([0.05 0.95], length(ec90_2(:,2))-1);
ts90_3 = tinv([0.05 0.95], length(ec90_3(:,2))-1);


%Confidence Interval
CI90_1 = mean(ec90_1(:,2)) + ts90_1*SEM90_1;
CI90_2 = mean(ec90_2(:,2)) + ts90_2*SEM90_2;
CI90_3 = mean(ec90_3(:,2)) + ts90_3*SEM90_3;

%Standard deviation
std90_1 = std(ec90_1(:,2));
std90_2 = std(ec90_2(:,2));
std90_3 = std(ec90_3(:,2));


%Interquantile range
r90_1 = iqr(ec90_1(:,2));
r90_2 = iqr(ec90_2(:,2));
r90_3 = iqr(ec90_3(:,2));


%Mean
mean_ec90_1 = mean(ec90_1(:,2));
mean_ec90_2 = mean(ec90_2(:,2));
mean_ec90_3 = mean(ec90_3(:,2));

%Variance
var_ec90_1 = var(ec90_1(:,2));
var_ec90_2 = var(ec90_2(:,2));
var_ec90_3 = var(ec90_3(:,2));


%ANOVA
all_ec90s = [ec90_1(:,2) ec90_2(:,2) ec90_3(:,2)];
group = ["rep1" "rep2" "rep3"];
[p_90,tbl_90,stats_90] = anova1(all_ec90s, group);
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
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_alltable', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_all', '.jpg']);
close all

ec90_1_2 = [ec90_1(:,2) ec90_2(:,2)];
group = ["rep1" "rep2"];
[p_90_1_2, tbl_90_1_2, stats_90_1_2] = anova1(ec90_1_2, group);
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
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep2_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep2', '.jpg']);
close all

ec90_1_3 = [ec90_1(:,2) ec90_3(:,2)];
group = ["rep1" "rep3"];
[p_90_1_3, tbl_90_1_3, stats_90_1_3] = anova1(ec90_1_3, group);
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
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep1_rep3', '.jpg']);
close all

ec90_2_3 = [ec90_2(:,2) ec90_3(:,2)];
group = ["rep2" "rep3"];
[p_90_2_3, tbl_90_2_3, stats_90_2_3] = anova1(ec90_2_3, group);
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
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep2_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_EC90_rep2_rep3', '.jpg']);
close all


%Coefficient of Variation
CV90_1 = (std(ec90_1(:,2))/mean(ec90_1(:,2))) * 100;
CV90_2 = (std(ec90_2(:,2))/mean(ec90_2(:,2))) * 100;
CV90_3 = (std(ec90_3(:,2))/mean(ec90_3(:,2))) * 100;


%Histograms
subplot(1,3,1)
hist_ec90_1 = histogram(ec90_1(:,2));
title('EC50 repression coefficient 1 - randomly drawn parameters')
subplot(1,3,2)
hist_ec90_2 = histogram(ec90_2(:,2));
title('EC50 repression coefficient 2 - randomly drawn parameters')
subplot(1,3,3)
hist_ec90_3 = histogram(ec90_3(:,2));
title('EC50 repression coefficient 3 - randomly drawn parameters')




%Check if the order of dose response curves for the same parameters but
%with different repression coefficient is the same all the time.


%save in one structure

statistics_TetRTup1_EC90.sem1 = SEM90_1;
statistics_TetRTup1_EC90.sem2 = SEM90_2;
statistics_TetRTup1_EC90.sem3 = SEM90_3;

statistics_TetRTup1_EC90.tscore1 = ts90_1;
statistics_TetRTup1_EC90.tscore2 = ts90_2;
statistics_TetRTup1_EC90.tscore3 = ts90_3;

statistics_TetRTup1_EC90.CI1 = CI90_1;
statistics_TetRTup1_EC90.CI2 = CI90_2;
statistics_TetRTup1_EC90.CI3 = CI90_3;

statistics_TetRTup1_EC90.std1 = std90_1;
statistics_TetRTup1_EC90.std2 = std90_2;
statistics_TetRTup1_EC90.std3 = std90_3;

statistics_TetRTup1_EC90.interquartile_range1 = r90_1;
statistics_TetRTup1_EC90.interquartile_range2 = r90_2;
statistics_TetRTup1_EC90.interquartile_range3 = r90_3;

statistics_TetRTup1_EC90.mean1 = mean_ec90_1;
statistics_TetRTup1_EC90.mean2 = mean_ec90_2;
statistics_TetRTup1_EC90.mean3 = mean_ec90_3;

statistics_TetRTup1_EC90.var1 = var_ec90_1;
statistics_TetRTup1_EC90.var2 = var_ec90_2;
statistics_TetRTup1_EC90.var3 = var_ec90_3;

statistics_TetRTup1_EC90.ANOVA.all.p = p_90;
statistics_TetRTup1_EC90.ANOVA.all.tbl = tbl_90;
statistics_TetRTup1_EC90.ANOVA.all.stats = stats_90;

statistics_TetRTup1_EC90.ANOVA.rep1_rep3.p = p_90_1_3;
statistics_TetRTup1_EC90.ANOVA.rep1_rep3.tbl = tbl_90_1_3;
statistics_TetRTup1_EC90.ANOVA.rep1_rep3.stats = stats_90_1_3;

statistics_TetRTup1_EC90.ANOVA.rep1_rep2.p = p_90_1_2;
statistics_TetRTup1_EC90.ANOVA.rep1_rep2.tbl = tbl_90_1_2;
statistics_TetRTup1_EC90.ANOVA.rep1_rep2.stats = stats_90_1_2;

statistics_TetRTup1_EC90.ANOVA.rep2_rep3.p = p_90_2_3;
statistics_TetRTup1_EC90.ANOVA.rep2_rep3.tbl = tbl_90_2_3;
statistics_TetRTup1_EC90.ANOVA.rep2_rep3.stats = stats_90_2_3;

statistics_TetRTup1_EC90.CV1 = CV90_1;
statistics_TetRTup1_EC90.CV2 = CV90_2;
statistics_TetRTup1_EC90.CV3 = CV90_3;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/',datestr(now, 'dd-mmm-yyyy'),'_statistics_TetRTup1_EC90.mat'],'statistics_TetRTup1_EC90')




