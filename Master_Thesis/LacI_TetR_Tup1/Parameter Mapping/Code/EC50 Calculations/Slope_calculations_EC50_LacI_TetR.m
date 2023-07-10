%Calculations of slope around specific points for the parameters drawn from
%hyperspace
load("27-Apr-2023ec50_1_TetRTup1.mat");
load("27-Apr-2023ec50_2_TetRTup1.mat");
load("27-Apr-2023ec50_3_TetRTup1.mat");
load("27-Apr-2023SimFluoValues1_LacI_TetRTup1_EC50.mat");
load("27-Apr-2023SimFluoValues2_LacI_TetRTup1_EC50.mat");
load("27-Apr-2023SimFluoValues3_LacI_TetRTup1_EC50.mat");


%Slopes around EC50
for k = 1:length(ec50_1)
    simfluovalue = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k);
    doses = log10(SimFluoValues1combined.Dose(SimFluoValues1combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec50_1(k,3)));
    lower_fluovalue = simfluovalue(closest_index - 1);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index-1);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec50_1(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC50_1_LacI_TetR_Tup1', '.mat'], 'slopes_ec50_1');

for k = 1:length(ec50_2)
    simfluovalue = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k);
    doses = log10(SimFluoValues2combined.Dose(SimFluoValues2combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec50_2(k,3)));
    lower_fluovalue = simfluovalue(closest_index - 1);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index-1);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec50_2(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC50_2_LacI_TetR_Tup1', '.mat'], 'slopes_ec50_2');

for k = 1:length(ec50_3)
    simfluovalue = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    doses = log10(SimFluoValues3combined.Dose(SimFluoValues3combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec50_3(k,3)));
    lower_fluovalue = simfluovalue(closest_index - 1);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index-1);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec50_3(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC50_3_LacI_TetR_Tup1', '.mat'], 'slopes_ec50_3');

%statistics for slopes around EC50
%standard error
SEM_slopes_50_1 = std(slopes_ec50_1(:))/sqrt(length(slopes_ec50_1(:)));
SEM_slopes_50_2 = std(slopes_ec50_2(:))/sqrt(length(slopes_ec50_2(:)));
SEM_slopes_50_3 = std(slopes_ec50_3(:))/sqrt(length(slopes_ec50_3(:)));

%T-Score
ts_slopes_50_1 = tinv([0.05 0.95], length(slopes_ec50_1(:))-1);
ts_slopes_50_2 = tinv([0.05 0.95], length(slopes_ec50_2(:))-1);
ts_slopes_50_3 = tinv([0.05 0.95], length(slopes_ec50_3(:))-1);

%Confidence Interval
CI_slopes_50_1 = mean(slopes_ec50_1(:)) + ts_slopes_50_1*SEM_slopes_50_1;
CI_slopes_50_2 = mean(slopes_ec50_2(:)) + ts_slopes_50_2*SEM_slopes_50_2;
CI_slopes_50_3 = mean(slopes_ec50_3(:)) + ts_slopes_50_3*SEM_slopes_50_3;

%Standard deviation
std_slopes_50_1 = std(slopes_ec50_1(:));
std_slopes_50_2 = std(slopes_ec50_2(:));
std_slopes_50_3 = std(slopes_ec50_3(:));

%Interquantile range
r_slopes_50_1 = iqr(slopes_ec50_1(:));
r_slopes_50_2 = iqr(slopes_ec50_2(:));
r_slopes_50_3 = iqr(slopes_ec50_3(:));

%Mean
mean_slopes_ec50_1 = mean(slopes_ec50_1(:));
mean_slopes_ec50_2 = mean(slopes_ec50_2(:));
mean_slopes_ec50_3 = mean(slopes_ec50_3(:));

%Variance
var_slopes_ec50_1 = var(slopes_ec50_1(:));
var_slopes_ec50_2 = var(slopes_ec50_2(:));
var_slopes_ec50_3 = var(slopes_ec50_3(:));

%ANOVA
all_slopes_ec50s = [slopes_ec50_1(:) slopes_ec50_2(:) slopes_ec50_3(:)];
group = ["rep1" "rep2" "rep3"];
[p_50_slopes,tbl_50_slopes,stats_50_slopes] = anova1(all_slopes_ec50s, group);
title('ANOVA slopes around EC50 - all Rep', 'FontSize',70);
ax = gca;
ax.FontSize = 40;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
ax = gca; % Get the current axes
ax.XAxis.FontSize = 40; % X axis font size
ax.YAxis.FontSize = 40; % Y axis font size
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_all_table_TetRTup1', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_all_TetRTup1', '.jpg']);
close all

ec50_slopes_1_2 = [slopes_ec50_1(:) slopes_ec50_2(:)];
group = ["rep1" "rep2"];
[p_50_1_2_slopes, tbl_50_1_2_slopes, stats_50_1_2_slopes] = anova1(ec50_slopes_1_2, group);
title('ANOVA slopes around EC50 - 1st and 2nd repression coefficient', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
ax = gca; % Get the current axes
ax.XAxis.FontSize = 40; % X axis font size
ax.YAxis.FontSize = 40; % Y axis font size
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep1_rep2_table_TetRTup1', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep1_rep2_TetRTup1', '.jpg']);
close all

ec50_slopes_1_3 = [slopes_ec50_1(:) slopes_ec50_3(:)];
group = ["rep1" "rep3"];
[p_50_1_3_slopes, tbl_50_1_3_slopes, stats_50_1_3_slopes] = anova1(ec50_slopes_1_3, group);
title('ANOVA slopes around EC50 - 1st and 3rd repression coefficient', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep1_rep3_table_TetRTup1', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep1_rep3_TetRTup1', '.jpg']);
close all

ec50_slopes_2_3 = [slopes_ec50_2(:) slopes_ec50_3(:)];
group = ["rep2" "rep3"];
[p_50_2_3_slopes, tbl_50_2_3_slopes, stats_50_2_3_slopes] = anova1(ec50_slopes_2_3, group);
title('ANOVA slopes around EC50 - 2nd and 3rd repression coefficient', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep2_rep3_table_TetRTup1', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep2_rep3_TetRTup1', '.jpg']);
close all

%Coefficient of Variation
CV50_1_slopes = (std(slopes_ec50_1(:))/mean(slopes_ec50_1(:))) * 100;
CV50_2_slopes = (std(slopes_ec50_2(:))/mean(slopes_ec50_2(:))) * 100;
CV50_3_slopes = (std(slopes_ec50_3(:))/mean(slopes_ec50_3(:))) * 100;

%Histograms
subplot(1,3,1)
hist_slopes_ec50_1 = histogram(slopes_ec50_1(:));
title('Slopes around EC50 repression coefficient 1 - LacI & TetRTup1')

subplot(1,3,2)
hist_slopes_ec50_2 = histogram(slopes_ec50_2(:));
title('Slopes around EC50 repression coefficient 2 - LacI & TetRTup1')

subplot(1,3,3)
hist_slopes_ec50_3 = histogram(slopes_ec50_3(:));
title('Slopes around EC50 repression coefficient 3 - LacI & TetRTup1')
set(figure(1), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/', datestr(now, 'dd-mmm-yyyy'),'histograms_slopes_EC50', '.jpg']);
close all



%save in one structure

statistics_LacI_TetRTup1_slopes.sem1 = SEM_slopes_50_1;
statistics_LacI_TetRTup1_slopes.sem2 = SEM_slopes_50_2;
statistics_LacI_TetRTup1_slopes.sem3 = SEM_slopes_50_3;

statistics_LacI_TetRTup1_slopes.tscore1 = ts_slopes_50_1;
statistics_LacI_TetRTup1_slopes.tscore2 = ts_slopes_50_2;
statistics_LacI_TetRTup1_slopes.tscore3 = ts_slopes_50_3;

statistics_LacI_TetRTup1_slopes.CI1 = CI_slopes_50_1;
statistics_LacI_TetRTup1_slopes.CI2 = CI_slopes_50_2;
statistics_LacI_TetRTup1_slopes.CI3 = CI_slopes_50_3;

statistics_LacI_TetRTup1_slopes.std1 = std_slopes_50_1;
statistics_LacI_TetRTup1_slopes.std2 = std_slopes_50_2;
statistics_LacI_TetRTup1_slopes.std3 = std_slopes_50_3;

statistics_LacI_TetRTup1_slopes.interquartile_range1 = r_slopes_50_1;
statistics_LacI_TetRTup1_slopes.interquartile_range2 = r_slopes_50_2;
statistics_LacI_TetRTup1_slopes.interquartile_range3 = r_slopes_50_3;

statistics_LacI_TetRTup1_slopes.mean1 = mean_slopes_ec50_1;
statistics_LacI_TetRTup1_slopes.mean2 = mean_slopes_ec50_2;
statistics_LacI_TetRTup1_slopes.mean3 = mean_slopes_ec50_3;

statistics_LacI_TetRTup1_slopes.var1 = var_slopes_ec50_1;
statistics_LacI_TetRTup1_slopes.var2 = var_slopes_ec50_2;
statistics_LacI_TetRTup1_slopes.var3 = var_slopes_ec50_3;

statistics_LacI_TetRTup1_slopes.ANOVA.all.p = p_50_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.all.tbl = tbl_50_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.all.stats = stats_50_slopes;

statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep3.p = p_50_1_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep3.tbl = tbl_50_1_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep3.stats = stats_50_1_3_slopes;

statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep2.p = p_50_1_2_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep2.tbl = tbl_50_1_2_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep2.stats = stats_50_1_2_slopes;

statistics_LacI_TetRTup1_slopes.ANOVA.rep2_rep3.p = p_50_2_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep2_rep3.tbl = tbl_50_2_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep2_rep3.stats = stats_50_2_3_slopes;

statistics_LacI_TetRTup1_slopes.CV1 = CV50_1_slopes;
statistics_LacI_TetRTup1_slopes.CV2 = CV50_2_slopes;
statistics_LacI_TetRTup1_slopes.CV3 = CV50_3_slopes;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations/',datestr(now, 'dd-mmm-yyyy'),'_statistics_LacI_TetRTup1_slopes.mat'],'statistics_LacI_TetRTup1_slopes')


