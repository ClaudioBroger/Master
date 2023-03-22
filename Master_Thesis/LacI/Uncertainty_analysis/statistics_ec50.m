%standard error
SEM1 = std(ec50_1(:,2))/length(ec50_1(:,2));
SEM2 = std(ec50_2(:,2))/length(ec50_2(:,2));
SEM3 = std(ec50_3(:,2))/length(ec50_3(:,2));

%T-Score
ts1 = tinv([0.05 0.95], length(ec50_1(:,2))-1);
ts2 = tinv([0.05 0.95], length(ec50_2(:,2))-1);
ts3 = tinv([0.05 0.95], length(ec50_3(:,2))-1);

%Confidence Interval
CI1 = mean(ec50_1(:,2)) + ts1*SEM1;
CI2 = mean(ec50_2(:,2)) + ts2*SEM2;
CI3 = mean(ec50_3(:,2)) + ts3*SEM3;

sorted_ec50 = sort(ec50_1(:,2));

%ANOVA
mean_ec50_1 = mean(ec50_1(:,2));
mean_ec50_2 = mean(ec50_2(:,2));
mean_ec50_3 = mean(ec50_3(:,2));

var_ec50_1 = var(ec50_1(:,2));
var_ec50_2 = var(ec50_2(:,2));
var_ec50_3 = var(ec50_3(:,2));

mean_all = [mean_ec50_1 mean_ec50_2 mean_ec50_3];
group = ["rep1" "rep2" "rep3"];

[p_mean,tbl_mean,stats_mean] = anova1(mean_all,group);

var_all = [var_ec50_1 var_ec50_2 var_ec50_3];
[p_var,tbl_var,stats_var] = anova1(var_all,group);

%Coefficient of Variation
CV1 = (std(ec50_1(:,2))/mean(ec50_1(:,2))) * 100;
CV2 = (std(ec50_2(:,2))/mean(ec50_2(:,2))) * 100;
CV3 = (std(ec50_3(:,2))/mean(ec50_3(:,2))) * 100;

