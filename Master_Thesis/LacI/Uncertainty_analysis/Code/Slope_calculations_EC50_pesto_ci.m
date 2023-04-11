%Calculations of slope around specific points for the parameters drawn from
%PESTO CI

%Slopes around EC50
for k = 1:length(ec50_1_pesto)
    simfluovalue = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k);
    doses = log10(SimFluoValues1combined.dose(SimFluoValues1combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec50_1_pesto(k,3)));
    lower_fluovalue = simfluovalue(closest_index - 1);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index-1);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec50_1_pesto(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC50_1_pesto_ci', '.mat'], 'slopes_ec50_1_pesto');

for k = 1:length(ec50_2_pesto)
    simfluovalue = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k);
    doses = log10(SimFluoValues2combined.dose(SimFluoValues2combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec50_2_pesto(k,3)));
    lower_fluovalue = simfluovalue(closest_index - 1);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index-1);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec50_2_pesto(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC50_2_pesto_ci', '.mat'], 'slopes_ec50_2_pesto');

for k = 1:length(ec50_3_pesto)
    simfluovalue = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    doses = log10(SimFluoValues3combined.dose(SimFluoValues3combined.Draw == k));
    [minvalue closest_index] = min(abs(simfluovalue - ec50_3_pesto(k,3)));
    lower_fluovalue = simfluovalue(closest_index - 1);
    upper_fluovalue = simfluovalue(closest_index + 1);
    deltafluo = upper_fluovalue - lower_fluovalue;
    lower_dose = doses(closest_index-1);
    upper_dose = doses(closest_index+1);
    deltadose = upper_dose - lower_dose;
    slopes_ec50_3_pesto(k,:) = deltafluo/deltadose;
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'slopes_EC50_3_pesto_ci', '.mat'], 'slopes_ec50_3_pesto');

%statistics for slopes around EC50
%standard error
SEM_slopes_50_1_pesto = std(slopes_ec50_1_pesto(:))/sqrt(length(slopes_ec50_1_pesto(:)));
SEM_slopes_50_2_pesto = std(slopes_ec50_2_pesto(:))/sqrt(length(slopes_ec50_2_pesto(:)));
SEM_slopes_50_3_pesto = std(slopes_ec50_3_pesto(:))/sqrt(length(slopes_ec50_3_pesto(:)));

%T-Score
ts_slopes_50_1_pesto = tinv([0.05 0.95], length(slopes_ec50_1_pesto(:))-1);
ts_slopes_50_2_pesto = tinv([0.05 0.95], length(slopes_ec50_2_pesto(:))-1);
ts_slopes_50_3_pesto = tinv([0.05 0.95], length(slopes_ec50_3_pesto(:))-1);

%Confidence Interval
CI_slopes_50_1_pesto = mean(slopes_ec50_1_pesto(:)) + ts_slopes_50_1_pesto*SEM_slopes_50_1_pesto;
CI_slopes_50_2_pesto = mean(slopes_ec50_2_pesto(:)) + ts_slopes_50_2_pesto*SEM_slopes_50_2_pesto;
CI_slopes_50_3_pesto = mean(slopes_ec50_3_pesto(:)) + ts_slopes_50_3_pesto*SEM_slopes_50_3_pesto;

%Standard deviation
std_slopes_50_1_pesto = std(slopes_ec50_1_pesto(:));
std_slopes_50_2_pesto = std(slopes_ec50_2_pesto(:));
std_slopes_50_3_pesto = std(slopes_ec50_3_pesto(:));

%Interquantile range
r_slopes_50_1_pesto = iqr(slopes_ec50_1_pesto(:));
r_slopes_50_2_pesto = iqr(slopes_ec50_2_pesto(:));
r_slopes_50_3_pesto = iqr(slopes_ec50_3_pesto(:));

%Mean
mean_slopes_ec50_1_pesto = mean(slopes_ec50_1_pesto(:));
mean_slopes_ec50_2_pesto = mean(slopes_ec50_2_pesto(:));
mean_slopes_ec50_3_pesto = mean(slopes_ec50_3_pesto(:));

%Variance
var_slopes_ec50_1_pesto = var(slopes_ec50_1_pesto(:));
var_slopes_ec50_2_pesto = var(slopes_ec50_2_pesto(:));
var_slopes_ec50_3_pesto = var(slopes_ec50_3_pesto(:));

%ANOVA
all_slopes_ec50s_pesto = [slopes_ec50_1_pesto(:) slopes_ec50_2_pesto(:) slopes_ec50_3_pesto(:)];
group = ["rep1" "rep2" "rep3"];
[p_50_slopes_pesto,tbl_50_slopes_pesto,stats_50_slopes_pesto] = anova1(all_slopes_ec50s_pesto, group);
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_all_table_pesto_ci', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_all_pesto_ci', '.jpg']);
close all

ec50_slopes_1_2_pesto = [slopes_ec50_1_pesto(:) slopes_ec50_2_pesto(:)];
group = ["rep1" "rep2"];
[p_50_1_2_slopes_pesto, tbl_50_1_2_slopes_pesto, stats_50_1_2_slopes_pesto] = anova1(ec50_slopes_1_2_pesto, group);
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep1_rep2_table_pesto_ci', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep1_rep2_pesto_ci', '.jpg']);
close all

ec50_slopes_1_3_pesto = [slopes_ec50_1_pesto(:) slopes_ec50_3_pesto(:)];
group = ["rep1" "rep3"];
[p_50_1_3_slopes_pesto, tbl_50_1_3_slopes_pesto, stats_50_1_3_slopes_pesto] = anova1(ec50_slopes_1_3_pesto, group);
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep1_rep3_table_pesto_ci', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep1_rep3_pesto_ci', '.jpg']);
close all

ec50_slopes_2_3_pesto = [slopes_ec50_2_pesto(:) slopes_ec50_3_pesto(:)];
group = ["rep2" "rep3"];
[p_50_2_3_slopes_pesto, tbl_50_2_3_slopes_pesto, stats_50_2_3_slopes_pesto] = anova1(ec50_slopes_2_3_pesto, group);
set(figure(2), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep2_rep3_table_pesto_ci', '.jpg']);
saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'ANOVA_slopes_EC50_rep2_rep3_pesto_ci', '.jpg']);
close all

%Coefficient of Variation
CV50_1_slopes_pesto = (std(slopes_ec50_1_pesto(:))/mean(slopes_ec50_1_pesto(:))) * 100;
CV50_2_slopes_pesto = (std(slopes_ec50_2_pesto(:))/mean(slopes_ec50_2_pesto(:))) * 100;
CV50_3_slopes_pesto = (std(slopes_ec50_3_pesto(:))/mean(slopes_ec50_3_pesto(:))) * 100;

%Histograms
subplot(1,3,1)
hist_slopes_ec50_1_pesto = histogram(slopes_ec50_1_pesto(:));
title('Slopes around EC50 repression coefficient 1 - parameters from PESTO')

subplot(1,3,2)
hist_slopes_ec50_2_pesto = histogram(slopes_ec50_2_pesto(:));
title('Slopes around EC50 repression coefficient 2 - parameters from PESTO')

subplot(1,3,3)
hist_slopes_ec50_3_pesto = histogram(slopes_ec50_3_pesto(:));
title('Slopes around EC50 repression coefficient 3 - parameters from PESTO')
set(figure(1), 'Position', get(0, 'Screensize'))
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/', datestr(now, 'dd-mmm-yyyy'),'histograms_slopes_EC50_pesto_ci', '.jpg']);
close all



%save in one structure

statistics_pesto_slopes.sem1 = SEM_slopes_50_1_pesto;
statistics_pesto_slopes.sem2 = SEM_slopes_50_2_pesto;
statistics_pesto_slopes.sem3 = SEM_slopes_50_3_pesto;

statistics_pesto_slopes.tscore1 = ts_slopes_50_1_pesto;
statistics_pesto_slopes.tscore2 = ts_slopes_50_2_pesto;
statistics_pesto_slopes.tscore3 = ts_slopes_50_3_pesto;

statistics_pesto_slopes.CI1 = CI_slopes_50_1_pesto;
statistics_pesto_slopes.CI2 = CI_slopes_50_2_pesto;
statistics_pesto_slopes.CI3 = CI_slopes_50_3_pesto;

statistics_pesto_slopes.std1 = std_slopes_50_1_pesto;
statistics_pesto_slopes.std2 = std_slopes_50_2_pesto;
statistics_pesto_slopes.std3 = std_slopes_50_3_pesto;

statistics_pesto_slopes.interquartile_range1 = r_slopes_50_1_pesto;
statistics_pesto_slopes.interquartile_range2 = r_slopes_50_2_pesto;
statistics_pesto_slopes.interquartile_range3 = r_slopes_50_3_pesto;

statistics_pesto_slopes.mean1 = mean_slopes_ec50_1_pesto;
statistics_pesto_slopes.mean2 = mean_slopes_ec50_2_pesto;
statistics_pesto_slopes.mean3 = mean_slopes_ec50_3_pesto;

statistics_pesto_slopes.var1 = var_slopes_ec50_1_pesto;
statistics_pesto_slopes.var2 = var_slopes_ec50_2_pesto;
statistics_pesto_slopes.var3 = var_slopes_ec50_3_pesto;

statistics_pesto_slopes.ANOVA.all.p = p_50_slopes_pesto;
statistics_pesto_slopes.ANOVA.all.tbl = tbl_50_slopes_pesto;
statistics_pesto_slopes.ANOVA.all.stats = stats_50_slopes_pesto;

statistics_pesto_slopes.ANOVA.rep1_rep3.p = p_50_1_3_slopes_pesto;
statistics_pesto_slopes.ANOVA.rep1_rep3.tbl = tbl_50_1_3_slopes_pesto;
statistics_pesto_slopes.ANOVA.rep1_rep3.stats = stats_50_1_3_slopes_pesto;

statistics_pesto_slopes.ANOVA.rep1_rep2.p = p_50_1_2_slopes_pesto;
statistics_pesto_slopes.ANOVA.rep1_rep2.tbl = tbl_50_1_2_slopes_pesto;
statistics_pesto_slopes.ANOVA.rep1_rep2.stats = stats_50_1_2_slopes_pesto;

statistics_pesto_slopes.ANOVA.rep2_rep3.p = p_50_2_3_slopes_pesto;
statistics_pesto_slopes.ANOVA.rep2_rep3.tbl = tbl_50_2_3_slopes_pesto;
statistics_pesto_slopes.ANOVA.rep2_rep3.stats = stats_50_2_3_slopes_pesto;

statistics_pesto_slopes.CV1 = CV50_1_slopes_pesto;
statistics_pesto_slopes.CV2 = CV50_2_slopes_pesto;
statistics_pesto_slopes.CV3 = CV50_3_slopes_pesto;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/',datestr(now, 'dd-mmm-yyyy'),'_statistics_pesto_ci_slopes.mat'],'statistics_pesto_slopes')

