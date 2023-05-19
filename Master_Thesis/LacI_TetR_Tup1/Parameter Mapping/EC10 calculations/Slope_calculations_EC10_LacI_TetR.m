%Calculations of slope around specific points for the parameters drawn from
%hyperspace
load("28-Apr-2023ec10_1_TetRTup1.mat");
load("28-Apr-2023ec10_2_TetRTup1.mat");
load("29-Apr-2023ec10_3_TetRTup1.mat");
load("27-Apr-2023SimFluoValues1_LacI_TetRTup1_EC50.mat");
load("27-Apr-2023SimFluoValues2_LacI_TetRTup1_EC50.mat");
load("27-Apr-2023SimFluoValues3_LacI_TetRTup1_EC50.mat");


%Slopes around EC50
for k = 1:length(ec10_1)
    simfluovalue = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k);
    doses = log10(SimFluoValues1combined.Dose(SimFluoValues1combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec10_1(k,3)));
    lower_fluovalue = simfluovalue(closest_index - 1);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index-1);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec10_1(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC10_1_LacI_TetR_Tup1', '.mat'], 'slopes_ec10_1');

for k = 1:length(ec10_2)
    simfluovalue = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k);
    doses = log10(SimFluoValues2combined.Dose(SimFluoValues2combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec10_2(k,3)));
    lower_fluovalue = simfluovalue(closest_index - 1);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index-1);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec10_2(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC10_2_LacI_TetR_Tup1', '.mat'], 'slopes_ec10_2');

for k = 1:length(ec10_3)
    simfluovalue = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    doses = log10(SimFluoValues3combined.Dose(SimFluoValues3combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec10_3(k,3)));
    lower_fluovalue = simfluovalue(closest_index - 1);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index-1);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec10_3(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC10_3_LacI_TetR_Tup1', '.mat'], 'slopes_ec10_3');

%statistics for slopes around EC50
%standard error
SEM_slopes_10_1 = std(slopes_ec10_1(:))/sqrt(length(slopes_ec10_1(:)));
SEM_slopes_10_2 = std(slopes_ec10_2(:))/sqrt(length(slopes_ec10_2(:)));
SEM_slopes_10_3 = std(slopes_ec10_3(:))/sqrt(length(slopes_ec10_3(:)));

%T-Score
ts_slopes_10_1 = tinv([0.05 0.95], length(slopes_ec10_1(:))-1);
ts_slopes_10_2 = tinv([0.05 0.95], length(slopes_ec10_2(:))-1);
ts_slopes_10_3 = tinv([0.05 0.95], length(slopes_ec10_3(:))-1);

%Confidence Interval
CI_slopes_10_1 = mean(slopes_ec10_1(:)) + ts_slopes_10_1*SEM_slopes_10_1;
CI_slopes_10_2 = mean(slopes_ec10_2(:)) + ts_slopes_10_2*SEM_slopes_10_2;
CI_slopes_10_3 = mean(slopes_ec10_3(:)) + ts_slopes_10_3*SEM_slopes_10_3;

%Standard deviation
std_slopes_10_1 = std(slopes_ec10_1(:));
std_slopes_10_2 = std(slopes_ec10_2(:));
std_slopes_10_3 = std(slopes_ec10_3(:));

%Interquantile range
r_slopes_10_1 = iqr(slopes_ec10_1(:));
r_slopes_10_2 = iqr(slopes_ec10_2(:));
r_slopes_10_3 = iqr(slopes_ec10_3(:));

%Mean
mean_slopes_ec10_1 = mean(slopes_ec10_1(:));
mean_slopes_ec10_2 = mean(slopes_ec10_2(:));
mean_slopes_ec10_3 = mean(slopes_ec10_3(:));

%Variance
var_slopes_ec10_1 = var(slopes_ec10_1(:));
var_slopes_ec10_2 = var(slopes_ec10_2(:));
var_slopes_ec10_3 = var(slopes_ec10_3(:));

%ANOVA
all_slopes_ec10s = [slopes_ec10_1(:) slopes_ec10_2(:) slopes_ec10_3(:)];
group = ["rep1" "rep2" "rep3"];
[p_10_slopes,tbl_10_slopes,stats_10_slopes] = anova1(all_slopes_ec10s, group);
figure(2)
hold on
title('ANOVA slopes around EC10 all repression coefficients - LacI & TetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC10_all_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC10_all', '.jpg']);
close all

ec10_slopes_1_2 = [slopes_ec10_1(:) slopes_ec10_2(:)];
group = ["rep1" "rep2"];
[p_10_1_2_slopes, tbl_10_1_2_slopes, stats_10_1_2_slopes] = anova1(ec10_slopes_1_2, group);
figure(2)
hold on
title('ANOVA slopes around EC10 repression coefficients 1 and 2 - LacI & TetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC10_rep1_rep2_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC10_rep1_rep2', '.jpg']);
close all

ec10_slopes_1_3 = [slopes_ec10_1(:) slopes_ec10_3(:)];
group = ["rep1" "rep3"];
[p_10_1_3_slopes, tbl_10_1_3_slopes, stats_10_1_3_slopes] = anova1(ec10_slopes_1_3, group);
figure(2)
hold on
title('ANOVA slopes around EC10 repression coefficients 1 and 3 - LacI & TetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC10_rep1_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC10_rep1_rep3', '.jpg']);
close all

ec10_slopes_2_3 = [slopes_ec10_2(:) slopes_ec10_3(:)];
group = ["rep2" "rep3"];
[p_10_2_3_slopes, tbl_10_2_3_slopes, stats_10_2_3_slopes] = anova1(ec10_slopes_2_3, group);
figure(2)
hold on
title('ANOVA slopes around EC10 repression coefficients 2 and 3 - LacI & TetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 20;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.5;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC10_rep2_rep3_table', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC10_rep2_rep3', '.jpg']);
close all

%Coefficient of Variation
CV10_1_slopes = (std(slopes_ec10_1(:))/mean(slopes_ec10_1(:))) * 100;
CV10_2_slopes = (std(slopes_ec10_2(:))/mean(slopes_ec10_2(:))) * 100;
CV10_3_slopes = (std(slopes_ec10_3(:))/mean(slopes_ec10_3(:))) * 100;

%Histograms
subplot(1,3,1)
hist_slopes_ec10_1 = histogram(slopes_ec10_1(:));
title('Slopes around EC10 repression coefficient 1 - LacI & TetRTup1')

subplot(1,3,2)
hist_slopes_ec10_2 = histogram(slopes_ec10_2(:));
title('Slopes around EC10 repression coefficient 2 - LacI & TetRTup1')

subplot(1,3,3)
hist_slopes_ec10_3 = histogram(slopes_ec10_3(:));
title('Slopes around EC10 repression coefficient 3 - LacI & TetRTup1')
set(figure(1), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'histograms_slopes_EC10', '.jpg']);
close all



%save in one structure

statistics_LacI_TetRTup1_slopes.sem1 = SEM_slopes_10_1;
statistics_LacI_TetRTup1_slopes.sem2 = SEM_slopes_10_2;
statistics_LacI_TetRTup1_slopes.sem3 = SEM_slopes_10_3;

statistics_LacI_TetRTup1_slopes.tscore1 = ts_slopes_10_1;
statistics_LacI_TetRTup1_slopes.tscore2 = ts_slopes_10_2;
statistics_LacI_TetRTup1_slopes.tscore3 = ts_slopes_10_3;

statistics_LacI_TetRTup1_slopes.CI1 = CI_slopes_10_1;
statistics_LacI_TetRTup1_slopes.CI2 = CI_slopes_10_2;
statistics_LacI_TetRTup1_slopes.CI3 = CI_slopes_10_3;

statistics_LacI_TetRTup1_slopes.std1 = std_slopes_10_1;
statistics_LacI_TetRTup1_slopes.std2 = std_slopes_10_2;
statistics_LacI_TetRTup1_slopes.std3 = std_slopes_10_3;

statistics_LacI_TetRTup1_slopes.interquartile_range1 = r_slopes_10_1;
statistics_LacI_TetRTup1_slopes.interquartile_range2 = r_slopes_10_2;
statistics_LacI_TetRTup1_slopes.interquartile_range3 = r_slopes_10_3;

statistics_LacI_TetRTup1_slopes.mean1 = mean_slopes_ec10_1;
statistics_LacI_TetRTup1_slopes.mean2 = mean_slopes_ec10_2;
statistics_LacI_TetRTup1_slopes.mean3 = mean_slopes_ec10_3;

statistics_LacI_TetRTup1_slopes.var1 = var_slopes_ec10_1;
statistics_LacI_TetRTup1_slopes.var2 = var_slopes_ec10_2;
statistics_LacI_TetRTup1_slopes.var3 = var_slopes_ec10_3;

statistics_LacI_TetRTup1_slopes.ANOVA.all.p = p_10_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.all.tbl = tbl_10_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.all.stats = stats_10_slopes;

statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep3.p = p_10_1_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep3.tbl = tbl_10_1_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep3.stats = stats_10_1_3_slopes;

statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep2.p = p_10_1_2_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep2.tbl = tbl_10_1_2_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep2.stats = stats_10_1_2_slopes;

statistics_LacI_TetRTup1_slopes.ANOVA.rep2_rep3.p = p_10_2_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep2_rep3.tbl = tbl_10_2_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep2_rep3.stats = stats_10_2_3_slopes;

statistics_LacI_TetRTup1_slopes.CV1 = CV10_1_slopes;
statistics_LacI_TetRTup1_slopes.CV2 = CV10_2_slopes;
statistics_LacI_TetRTup1_slopes.CV3 = CV10_3_slopes;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/',datestr(now, 'dd-mmm-yyyy'),'_statistics_LacI_TetRTup1_slopes.mat'],'statistics_LacI_TetRTup1_slopes')


