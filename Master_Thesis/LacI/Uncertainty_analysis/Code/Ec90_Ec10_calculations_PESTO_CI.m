%%Statistics EC90 PESTO CI

%SEM
SEM90_1_pesto_ci = std(ec90_1_pesto_ci(:,2))/sqrt(length(ec90_1_pesto_ci(:,2)));
SEM90_2_pesto_ci = std(ec90_2_pesto_ci(:,2))/sqrt(length(ec90_2_pesto_ci(:,2)));
SEM90_3_pesto_ci = std(ec90_3_pesto_ci(:,2))/sqrt(length(ec90_3_pesto_ci(:,2)));

%T-Score
ts90_1_pesto_ci = tinv([0.05 0.95], length(ec90_1_pesto_ci(:,2))-1);
ts90_2_pesto_ci = tinv([0.05 0.95], length(ec90_2_pesto_ci(:,2))-1);
ts90_3_pesto_ci = tinv([0.05 0.95], length(ec90_3_pesto_ci(:,2))-1);

%CI
CI90_1_pesto_ci = mean(ec90_1_pesto_ci(:,2)) + ts90_1_pesto_ci*SEM90_1_pesto_ci;
CI90_2_pesto_ci = mean(ec90_2_pesto_ci(:,2)) + ts90_2_pesto_ci*SEM90_2_pesto_ci;
CI90_3_pesto_ci = mean(ec90_3_pesto_ci(:,2)) + ts90_3_pesto_ci*SEM90_3_pesto_ci;

%STD
std90_1_pesto_ci = std(ec90_1_pesto_ci(:,2));
std90_2_pesto_ci = std(ec90_2_pesto_ci(:,2));
std90_3_pesto_ci = std(ec90_3_pesto_ci(:,2));

%interquantile range
r90_1_pesto_ci = iqr(ec90_1_pesto_ci(:,2));
r90_2_pesto_ci = iqr(ec90_2_pesto_ci(:,2));
r90_3_pesto_ci = iqr(ec90_3_pesto_ci(:,2));

%Mean
mean_ec90_1_pesto_ci = mean(ec90_1_pesto_ci(:,2));
mean_ec90_2_pesto_ci = mean(ec90_2_pesto_ci(:,2));
mean_ec90_3_pesto_ci = mean(ec90_3_pesto_ci(:,2));

%Var
var_ec90_1_pesto_ci = var(ec90_1_pesto_ci(:,2));
var_ec90_2_pesto_ci = var(ec90_2_pesto_ci(:,2));
var_ec90_3_pesto_ci = var(ec90_3_pesto_ci(:,2));

%ANOVA
all_ec90s_pesto_ci = [ec90_1_pesto_ci(:,2) ec90_2_pesto_ci(:,2) ec90_3_pesto_ci(:,2)];
group = ["rep1" "rep2" "rep3"];
[p_90_pesto_ci,tbl_90_pesto_ci,stats_90_pesto_ci] = anova1(all_ec90s_pesto_ci, group);

ec90_1_2_pesto_ci = [ec90_1_pesto_ci(:,2) ec90_2_pesto_ci(:,2)];
group = ["rep1" "rep2"];
[p_90_1_2_pesto_ci, tbl_90_1_2_pesto_ci, stats_90_1_2_pesto_ci] = anova1(ec90_1_2_pesto_ci, group);

ec90_1_3_pesto_ci = [ec90_1_pesto_ci(:,2) ec90_3_pesto_ci(:,2)];
group = ["rep1" "rep3"];
[p_90_1_3_pesto_ci, tbl_90_1_3_pesto_ci, stats_90_1_3_pesto_ci] = anova1(ec90_1_3_pesto_ci, group);

ec90_2_3_pesto_ci = [ec90_2_pesto_ci(:,2) ec90_3_pesto_ci(:,2)];
group = ["rep2" "rep3"];
[p_90_2_3_pesto_ci, tbl_90_2_3_pesto_ci, stats_90_2_3_pesto_ci] = anova1(ec90_2_3_pesto_ci, group);

%CV
CV90_1_pesto_ci = (std(ec90_1_pesto_ci(:,2))/mean(ec90_1_pesto_ci(:,2))) * 100;
CV90_2_pesto_ci = (std(ec90_2_pesto_ci(:,2))/mean(ec90_2_pesto_ci(:,2))) * 100;
CV90_3_pesto_ci = (std(ec90_3_pesto_ci(:,2))/mean(ec90_3_pesto_ci(:,2))) * 100;

%Histograms
subplot(1,3,1)
hist_ec90_1_pesto_ci = histogram(ec90_1_pesto_ci(:,2));
title('EC90 repression coefficient 1 - parameters from PESTO CI')
subplot(1,3,2)
hist_ec90_2_pesto_ci = histogram(ec90_2_pesto_ci(:,2));
title('EC90 repression coefficient 2 - parameters from PESTO Ci')
subplot(1,3,3)
hist_ec90_3_pesto_ci = histogram(ec90_3_pesto_ci(:,2));
title('EC90 repression coefficient 3 - parameters from PESTO Ci')

statistics_pesto_ci.sem1 = SEM90_1_pesto_ci;
statistics_pesto_ci.sem2 = SEM90_2_pesto_ci;
statistics_pesto_ci.sem3 = SEM90_3_pesto_ci;

statistics_pesto_ci.tscore1 = ts90_1_pesto_ci;
statistics_pesto_ci.tscore2 = ts90_2_pesto_ci;
statistics_pesto_ci.tscore3 = ts90_3_pesto_ci;

statistics_pesto_ci.CI1 = CI90_1_pesto_ci;
statistics_pesto_ci.CI2 = CI90_2_pesto_ci;
statistics_pesto_ci.CI3 = CI90_3_pesto_ci;

statistics_pesto_ci.std1 = std90_1_pesto_ci;
statistics_pesto_ci.std2 = std90_2_pesto_ci;
statistics_pesto_ci.std3 = std90_3_pesto_ci;

statistics_pesto_ci.interquartile_range1 = r90_1_pesto_ci;
statistics_pesto_ci.interquartile_range2 = r90_2_pesto_ci;
statistics_pesto_ci.interquartile_range3 = r90_3_pesto_ci;

statistics_pesto_ci.mean1 = mean_ec90_1_pesto_ci;
statistics_pesto_ci.mean2 = mean_ec90_2_pesto_ci;
statistics_pesto_ci.mean3 = mean_ec90_3_pesto_ci;

statistics_pesto_ci.var1 = var_ec90_1_pesto_ci;
statistics_pesto_ci.var2 = var_ec90_2_pesto_ci;
statistics_pesto_ci.var3 = var_ec90_3_pesto_ci;

statistics_pesto_ci.ANOVA.all.p = p_90_pesto_ci;
statistics_pesto_ci.ANOVA.all.tbl = tbl_90_pesto_ci;
statistics_pesto_ci.ANOVA.all.stats = stats_90_pesto_ci;

statistics_pesto_ci.ANOVA.rep1_rep3.p = p_90_1_3_pesto_ci;
statistics_pesto_ci.ANOVA.rep1_rep3.tbl = tbl_90_1_3_pesto_ci;
statistics_pesto_ci.ANOVA.rep1_rep3.stats = stats_90_1_3_pesto_ci;

statistics_pesto_ci.ANOVA.rep1_rep2.p = p_90_1_2_pesto_ci;
statistics_pesto_ci.ANOVA.rep1_rep2.tbl = tbl_90_1_2_pesto_ci;
statistics_pesto_ci.ANOVA.rep1_rep2.stats = stats_90_1_2_pesto_ci;

statistics_pesto_ci.ANOVA.rep2_rep3.p = p_90_2_3_pesto_ci;
statistics_pesto_ci.ANOVA.rep2_rep3.tbl = tbl_90_2_3_pesto_ci;
statistics_pesto_ci.ANOVA.rep2_rep3.stats = stats_90_2_3_pesto_ci;

statistics_pesto_ci.CV1 = CV90_1_pesto_ci;
statistics_pesto_ci.CV2 = CV90_2_pesto_ci;
statistics_pesto_ci.CV3 = CV90_3_pesto_ci;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/',datestr(now, 'dd-mmm-yyyy'),'_statistics__EC90_pesto_ci.mat'],'statistics_pesto_ci')

%%Statistics EC10

%SEM
SEM10_1_pesto_ci = std(ec10_1_pesto_ci(:,2))/sqrt(length(ec10_1_pesto_ci(:,2)));
SEM10_2_pesto_ci = std(ec10_2_pesto_ci(:,2))/sqrt(length(ec10_2_pesto_ci(:,2)));
SEM10_3_pesto_ci = std(ec10_3_pesto_ci(:,2))/sqrt(length(ec10_3_pesto_ci(:,2)));

%T-Score
ts10_1_pesto_ci = tinv([0.05 0.95], length(ec10_1_pesto_ci(:,2))-1);
ts10_2_pesto_ci = tinv([0.05 0.95], length(ec10_2_pesto_ci(:,2))-1);
ts10_3_pesto_ci = tinv([0.05 0.95], length(ec10_3_pesto_ci(:,2))-1);

%CI
CI10_1_pesto_ci = mean(ec10_1_pesto_ci(:,2)) + ts10_1_pesto_ci*SEM10_1_pesto_ci;
CI10_2_pesto_ci = mean(ec10_2_pesto_ci(:,2)) + ts10_2_pesto_ci*SEM10_2_pesto_ci;
CI10_3_pesto_ci = mean(ec10_3_pesto_ci(:,2)) + ts10_3_pesto_ci*SEM10_3_pesto_ci;

%STD
std10_1_pesto_ci = std(ec10_1_pesto_ci(:,2));
std10_2_pesto_ci = std(ec10_2_pesto_ci(:,2));
std10_3_pesto_ci = std(ec10_3_pesto_ci(:,2));

%interquantile range
r10_1_pesto_ci = iqr(ec10_1_pesto_ci(:,2));
r10_2_pesto_ci = iqr(ec10_2_pesto_ci(:,2));
r10_3_pesto_ci = iqr(ec10_3_pesto_ci(:,2));

%Mean
mean_ec10_1_pesto_ci = mean(ec10_1_pesto_ci(:,2));
mean_ec10_2_pesto_ci = mean(ec10_2_pesto_ci(:,2));
mean_ec10_3_pesto_ci = mean(ec10_3_pesto_ci(:,2));

%Var
var_ec10_1_pesto_ci = var(ec10_1_pesto_ci(:,2));
var_ec10_2_pesto_ci = var(ec10_2_pesto_ci(:,2));
var_ec10_3_pesto_ci = var(ec10_3_pesto_ci(:,2));

%ANOVA
all_ec10s_pesto_ci = [ec10_1_pesto_ci(:,2) ec10_2_pesto_ci(:,2) ec10_3_pesto_ci(:,2)];
group = ["rep1" "rep2" "rep3"];
[p_10_pesto_ci,tbl_10_pesto_ci,stats_10_pesto_ci] = anova1(all_ec10s_pesto_ci, group);

ec10_1_2_pesto_ci = [ec10_1_pesto_ci(:,2) ec10_2_pesto_ci(:,2)];
group = ["rep1" "rep2"];
[p_10_1_2_pesto_ci, tbl_10_1_2_pesto_ci, stats_10_1_2_pesto_ci] = anova1(ec10_1_2_pesto_ci, group);

ec10_1_3_pesto_ci = [ec10_1_pesto_ci(:,2) ec10_3_pesto_ci(:,2)];
group = ["rep1" "rep3"];
[p_10_1_3_pesto_ci, tbl_10_1_3_pesto_ci, stats_10_1_3_pesto_ci] = anova1(ec10_1_3_pesto_ci, group);

ec10_2_3_pesto_ci = [ec10_2_pesto_ci(:,2) ec10_3_pesto_ci(:,2)];
group = ["rep2" "rep3"];
[p_10_2_3_pesto_ci, tbl_10_2_3_pesto_ci, stats_10_2_3_pesto_ci] = anova1(ec10_2_3_pesto_ci, group);

%CV
CV10_1_pesto_ci = (std(ec10_1_pesto_ci(:,2))/mean(ec10_1_pesto_ci(:,2))) * 100;
CV10_2_pesto_ci = (std(ec10_2_pesto_ci(:,2))/mean(ec10_2_pesto_ci(:,2))) * 100;
CV10_3_pesto_ci = (std(ec10_3_pesto_ci(:,2))/mean(ec10_3_pesto_ci(:,2))) * 100;

%Histograms
subplot(1,3,1)
hist_ec10_1_pesto_ci = histogram(ec10_1_pesto_ci(:,2));
title('EC10 repression coefficient 1 - parameters from PESTO CI')
subplot(1,3,2)
hist_ec10_2_pesto_ci = histogram(ec10_2_pesto_ci(:,2));
title('EC10 repression coefficient 2 - parameters from PESTO CI')
subplot(1,3,3)
hist_ec10_3_pesto_ci = histogram(ec10_3_pesto_ci(:,2));
title('EC10 repression coefficient 3 - parameters from PESTO CI')

statistics_pesto_ci.sem1 = SEM10_1_pesto_ci;
statistics_pesto_ci.sem2 = SEM10_2_pesto_ci;
statistics_pesto_ci.sem3 = SEM10_3_pesto_ci;

statistics_pesto_ci.tscore1 = ts10_1_pesto_ci;
statistics_pesto_ci.tscore2 = ts10_2_pesto_ci;
statistics_pesto_ci.tscore3 = ts10_3_pesto_ci;

statistics_pesto_ci.CI1 = CI10_1_pesto_ci;
statistics_pesto_ci.CI2 = CI10_2_pesto_ci;
statistics_pesto_ci.CI3 = CI10_3_pesto_ci;

statistics_pesto_ci.std1 = std10_1_pesto_ci;
statistics_pesto_ci.std2 = std10_2_pesto_ci;
statistics_pesto_ci.std3 = std10_3_pesto_ci;

statistics_pesto_ci.interquartile_range1 = r10_1_pesto_ci;
statistics_pesto_ci.interquartile_range2 = r10_2_pesto_ci;
statistics_pesto_ci.interquartile_range3 = r10_3_pesto_ci;

statistics_pesto_ci.mean1 = mean_ec10_1_pesto_ci;
statistics_pesto_ci.mean2 = mean_ec10_2_pesto_ci;
statistics_pesto_ci.mean3 = mean_ec10_3_pesto_ci;

statistics_pesto_ci.var1 = var_ec10_1_pesto_ci;
statistics_pesto_ci.var2 = var_ec10_2_pesto_ci;
statistics_pesto_ci.var3 = var_ec10_3_pesto_ci;

statistics_pesto_ci.ANOVA.all.p = p_10_pesto_ci;
statistics_pesto_ci.ANOVA.all.tbl = tbl_10_pesto_ci;
statistics_pesto_ci.ANOVA.all.stats = stats_10_pesto_ci;

statistics_pesto_ci.ANOVA.rep1_rep3.p = p_10_1_3_pesto_ci;
statistics_pesto_ci.ANOVA.rep1_rep3.tbl = tbl_10_1_3_pesto_ci;
statistics_pesto_ci.ANOVA.rep1_rep3.stats = stats_10_1_3_pesto_ci;

statistics_pesto_ci.ANOVA.rep1_rep2.p = p_10_1_2_pesto_ci;
statistics_pesto_ci.ANOVA.rep1_rep2.tbl = tbl_10_1_2_pesto_ci;
statistics_pesto_ci.ANOVA.rep1_rep2.stats = stats_10_1_2_pesto_ci;

statistics_pesto_ci.ANOVA.rep2_rep3.p = p_10_2_3_pesto_ci;
statistics_pesto_ci.ANOVA.rep2_rep3.tbl = tbl_10_2_3_pesto_ci;
statistics_pesto_ci.ANOVA.rep2_rep3.stats = stats_10_2_3_pesto_ci;

statistics_pesto_ci.CV1 = CV10_1_pesto_ci;
statistics_pesto_ci.CV2 = CV10_2_pesto_ci;
statistics_pesto_ci.CV3 = CV10_3_pesto_ci;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO CI/',datestr(now, 'dd-mmm-yyyy'),'_statistics__EC10_PESTO_CI.mat'],'statistics_pesto_ci')
