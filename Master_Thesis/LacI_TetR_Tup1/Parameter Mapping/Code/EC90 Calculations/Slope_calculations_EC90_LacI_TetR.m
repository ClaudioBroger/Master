%Calculations of slope around specific points for the parameters drawn from
%hyperspace
load("28-Apr-2023ec90_1_TetRTup1.mat");
load("28-Apr-2023ec90_2_TetRTup1.mat");
load("28-Apr-2023ec90_3_TetRTup1.mat");
load("27-Apr-2023SimFluoValues1_LacI_TetRTup1_EC50.mat");
load("27-Apr-2023SimFluoValues2_LacI_TetRTup1_EC50.mat");
load("27-Apr-2023SimFluoValues3_LacI_TetRTup1_EC50.mat");


%Slopes around EC90
for k = 1:length(ec90_1)
    simfluovalue = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k);
    doses = log10(SimFluoValues1combined.Dose(SimFluoValues1combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec90_1(k,3)));
    lower_fluovalue = simfluovalue(closest_index);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec90_1(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC90_1_LacI_TetR_Tup1', '.mat'], 'slopes_ec90_1');

for k = 1:length(ec90_2)
    simfluovalue = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k);
    doses = log10(SimFluoValues2combined.Dose(SimFluoValues2combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec90_2(k,3)));
    lower_fluovalue = simfluovalue(closest_index);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec90_2(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC90_2_LacI_TetR_Tup1', '.mat'], 'slopes_ec90_2');

for k = 1:length(ec90_3)
    simfluovalue = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    doses = log10(SimFluoValues3combined.Dose(SimFluoValues3combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec90_3(k,3)));
    lower_fluovalue = simfluovalue(closest_index);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec90_3(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC90_3_LacI_TetR_Tup1', '.mat'], 'slopes_ec90_3');

%statistics for slopes around EC50
%standard error
SEM_slopes_90_1 = std(slopes_ec90_1(:))/sqrt(length(slopes_ec90_1(:)));
SEM_slopes_90_2 = std(slopes_ec90_2(:))/sqrt(length(slopes_ec90_2(:)));
SEM_slopes_90_3 = std(slopes_ec90_3(:))/sqrt(length(slopes_ec90_3(:)));

%T-Score
ts_slopes_90_1 = tinv([0.05 0.95], length(slopes_ec90_1(:))-1);
ts_slopes_90_2 = tinv([0.05 0.95], length(slopes_ec90_2(:))-1);
ts_slopes_90_3 = tinv([0.05 0.95], length(slopes_ec90_3(:))-1);

%Confidence Interval
CI_slopes_90_1 = mean(slopes_ec90_1(:)) + ts_slopes_90_1*SEM_slopes_90_1;
CI_slopes_90_2 = mean(slopes_ec90_2(:)) + ts_slopes_90_2*SEM_slopes_90_2;
CI_slopes_90_3 = mean(slopes_ec90_3(:)) + ts_slopes_90_3*SEM_slopes_90_3;

%Standard deviation
std_slopes_90_1 = std(slopes_ec90_1(:));
std_slopes_90_2 = std(slopes_ec90_2(:));
std_slopes_90_3 = std(slopes_ec90_3(:));

%Interquantile range
r_slopes_90_1 = iqr(slopes_ec90_1(:));
r_slopes_90_2 = iqr(slopes_ec90_2(:));
r_slopes_90_3 = iqr(slopes_ec90_3(:));

%Mean
mean_slopes_ec90_1 = mean(slopes_ec90_1(:));
mean_slopes_ec90_2 = mean(slopes_ec90_2(:));
mean_slopes_ec90_3 = mean(slopes_ec90_3(:));

%Variance
var_slopes_ec90_1 = var(slopes_ec90_1(:));
var_slopes_ec90_2 = var(slopes_ec90_2(:));
var_slopes_ec90_3 = var(slopes_ec90_3(:));

%ANOVA
all_slopes_ec90s = [slopes_ec90_1(:) slopes_ec90_2(:) slopes_ec90_3(:)];
group = ["rep1" "rep2" "rep3"];
[p_90_slopes,tbl_90_slopes,stats_90_slopes] = anova1(all_slopes_ec90s, group);
title('ANOVA slopes around EC90 - all Rep', 'FontSize',70);
ax = gca;
ax.FontSize = 40;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC90_all_table_TetRTup1', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC90_all_TetRTup1', '.jpg']);
close all

ec90_slopes_1_2 = [slopes_ec90_1(:) slopes_ec90_2(:)];
group = ["rep1" "rep2"];
[p_90_1_2_slopes, tbl_90_1_2_slopes, stats_90_1_2_slopes] = anova1(ec90_slopes_1_2, group);
title('ANOVA slopes around EC90 - 1st and 2nd repression coefficient - LacI & TetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC90_rep1_rep2_table_TetRTup1', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC90_rep1_rep2_TetRTup1', '.jpg']);
close all

ec90_slopes_1_3 = [slopes_ec90_1(:) slopes_ec90_3(:)];
group = ["rep1" "rep3"];
[p_90_1_3_slopes, tbl_90_1_3_slopes, stats_90_1_3_slopes] = anova1(ec90_slopes_1_3, group);
title('ANOVA slopes around EC90 - 1st and 3rd repression coefficient - LacI & TetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC90_rep1_rep3_table_TetRTup1', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC90_rep1_rep3_TetRTup1', '.jpg']);
close all

ec90_slopes_2_3 = [slopes_ec90_2(:) slopes_ec90_3(:)];
group = ["rep2" "rep3"];
[p_90_2_3_slopes, tbl_90_2_3_slopes, stats_90_2_3_slopes] = anova1(ec90_slopes_2_3, group);
title('ANOVA slopes around EC90 - 2nd and 3rd repression coefficient - LacI & TetRTup1', 'FontSize',24);
ax = gca;
ax.FontSize = 30;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
end
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC90_rep2_rep3_table_TetRTup1', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC90_rep2_rep3_TetRTup1', '.jpg']);
close all

%Coefficient of Variation
CV90_1_slopes = (std(slopes_ec90_1(:))/mean(slopes_ec90_1(:))) * 100;
CV90_2_slopes = (std(slopes_ec90_2(:))/mean(slopes_ec90_2(:))) * 100;
CV90_3_slopes = (std(slopes_ec90_3(:))/mean(slopes_ec90_3(:))) * 100;

%Histograms
subplot(1,3,1)
hist_slopes_ec90_1 = histogram(slopes_ec90_1(:));
title('Slopes around EC90 repression coefficient 1 - LacI & TetRTup1')

subplot(1,3,2)
hist_slopes_ec90_2 = histogram(slopes_ec90_2(:));
title('Slopes around EC90 repression coefficient 2 - LacI & TetRTup1')

subplot(1,3,3)
hist_slopes_ec90_3 = histogram(slopes_ec90_3(:));
title('Slopes around EC90 repression coefficient 3 - LacI & TetRTup1')
set(figure(1), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/', datestr(now, 'dd-mmm-yyyy'),'histograms_slopes_EC90', '.jpg']);
close all



%save in one structure

statistics_LacI_TetRTup1_slopes.sem1 = SEM_slopes_90_1;
statistics_LacI_TetRTup1_slopes.sem2 = SEM_slopes_90_2;
statistics_LacI_TetRTup1_slopes.sem3 = SEM_slopes_90_3;

statistics_LacI_TetRTup1_slopes.tscore1 = ts_slopes_90_1;
statistics_LacI_TetRTup1_slopes.tscore2 = ts_slopes_90_2;
statistics_LacI_TetRTup1_slopes.tscore3 = ts_slopes_90_3;

statistics_LacI_TetRTup1_slopes.CI1 = CI_slopes_90_1;
statistics_LacI_TetRTup1_slopes.CI2 = CI_slopes_90_2;
statistics_LacI_TetRTup1_slopes.CI3 = CI_slopes_90_3;

statistics_LacI_TetRTup1_slopes.std1 = std_slopes_90_1;
statistics_LacI_TetRTup1_slopes.std2 = std_slopes_90_2;
statistics_LacI_TetRTup1_slopes.std3 = std_slopes_90_3;

statistics_LacI_TetRTup1_slopes.interquartile_range1 = r_slopes_90_1;
statistics_LacI_TetRTup1_slopes.interquartile_range2 = r_slopes_90_2;
statistics_LacI_TetRTup1_slopes.interquartile_range3 = r_slopes_90_3;

statistics_LacI_TetRTup1_slopes.mean1 = mean_slopes_ec90_1;
statistics_LacI_TetRTup1_slopes.mean2 = mean_slopes_ec90_2;
statistics_LacI_TetRTup1_slopes.mean3 = mean_slopes_ec90_3;

statistics_LacI_TetRTup1_slopes.var1 = var_slopes_ec90_1;
statistics_LacI_TetRTup1_slopes.var2 = var_slopes_ec90_2;
statistics_LacI_TetRTup1_slopes.var3 = var_slopes_ec90_3;

statistics_LacI_TetRTup1_slopes.ANOVA.all.p = p_90_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.all.tbl = tbl_90_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.all.stats = stats_90_slopes;

statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep3.p = p_90_1_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep3.tbl = tbl_90_1_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep3.stats = stats_90_1_3_slopes;

statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep2.p = p_90_1_2_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep2.tbl = tbl_90_1_2_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep1_rep2.stats = stats_90_1_2_slopes;

statistics_LacI_TetRTup1_slopes.ANOVA.rep2_rep3.p = p_90_2_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep2_rep3.tbl = tbl_90_2_3_slopes;
statistics_LacI_TetRTup1_slopes.ANOVA.rep2_rep3.stats = stats_90_2_3_slopes;

statistics_LacI_TetRTup1_slopes.CV1 = CV90_1_slopes;
statistics_LacI_TetRTup1_slopes.CV2 = CV90_2_slopes;
statistics_LacI_TetRTup1_slopes.CV3 = CV90_3_slopes;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC90 calculations/',datestr(now, 'dd-mmm-yyyy'),'_statistics_LacI_TetRTup1_slopes.mat'],'statistics_LacI_TetRTup1_slopes')


