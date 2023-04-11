%PESTO uncertainty

load('0324_PestoResults_cleanup.mat')

dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');
%load model for uncertainty analysis
model = sbioloadproject('LacImodel_uncertainty');
num_draws = 200;
counter = 0;

[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));

parameters.name(1) = {"kLacI"};
parameters.name(5) = {"kCit"};
parameters.name(7) = {"dCit"};
parameters.name(8) = {"CitL"};
parameters.name(10) = {"LacIrep"};
parameters.name(11) = {"LacIrep2"};
parameters.name(12) = {"LacIrep3"};


parameters.name = cell2table(parameters.name);

[~,index] = ismember(ParaNames,parameters.name{:,:});
index = nonzeros(index);

min_par = parameters.min(index);
max_par = parameters.max(index);

CI = parameters.CI.S(index,:);

prob = parameters.S.par(index,:)';

for draw=1:num_draws
    for k=1:length(CI)
        pd = fitdist(prob(:,k),'Kernel');
         X = min_par(k):0.001:max_par(k);
         Y = pdf(pd,X);
         figure(1)
         subplot(2,6,k)
         plot(Y)
         rand_parameter(draw,k) = randsample(X,1,true,Y);
         
    end
end


% num_draws = 30;
ParaNames(6:7) = [];
% for k = 1:num_draws
%     for j = 1:length(CI)
%         rand_parameter(k,j) = CI(j,1) + (CI(j,2) - CI(j,1)).*rand;
%     end
% end

rand_parameter(:,1:6) = 10.^rand_parameter(:,1:6);
rand_parameter(:,8:11) = 10.^rand_parameter(:,8:11);

rand_parameter = array2table(rand_parameter);
rand_parameter.Properties.VariableNames = ParaNames;


rand_parameter.mu(:) = 0.0077;
rand_parameter.kmaturation(:) = 0.0173;


rand_parameter.nMperUnit(:) = 10^1.1;
[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
rand_parameter = rand_parameter(:,[string(ParaNames)]);
rand_parameter.dCit = zeros(num_draws,1);

%  save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/',datestr(now, 'dd-mmm-yyyy'),'_rand_parameter_uncertainty_analysis_ec50.mat'],'rand_parameter')

C = linspecer(num_draws);
for draw = 1:num_draws
    
    %define parameter values as the values randomly drawn
     paraValues = table2array(rand_parameter(draw,:));
     paraValues = paraValues';
    %TFlMode = FlMode;
    
     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;
    
        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;
            dataPos = strcat(dataPath, "data.mat");
            data1 = load(dataPos);
            data1.dose = [linspace(min(data1.dose), 10^3, 250) linspace(1+10^3, 10^7, 1500) linspace(1+10^7, max(data1.dose), 100)];
            data1.dose = data1.dose';
            
            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues1 = simulate_DR_IPTG(para,data1,ParaNames,model);
            DataMeans = data1.data.means;
            DataStd = data1.data.std;
            data1.dose(1) = 1;
        

        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;
            dataPos = strcat(dataPath, "data.mat");
            data2 = load(dataPos);
            data2.dose = [linspace(min(data1.dose), 10^5, 100) linspace(1+10^5, 10^7, 250) linspace(1+10^7, max(data1.dose), 100)];
            data2.dose = data2.dose';

            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep2', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues2 = simulate_DR_IPTG(para,data2,ParaNames,model);
            DataMeans = data2.data.means;
            DataStd = data2.data.std;
            data2.dose(1) = 1;

     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;
            dataPos = strcat(dataPath, "data.mat");
            data3 = load(dataPos);
            data3.dose = [linspace(min(data1.dose), 10^5, 100) linspace(1+10^5, 10^9, 250) linspace(1+10^9, 10^10, 100)];
            data3.dose = data3.dose';

            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues3 = simulate_DR_IPTG(para,data3,ParaNames,model);
            DataMeans = data3.data.means;
            DataStd = data3.data.std;
            data3.dose(1) = 1;

            %[hill3(draw,1) ec50_3(draw,1)] = doseResponse(data3.dose, SimFluoValues3);

            
        %set SimFluoValues to table
        SimFluoValues1 = array2table(SimFluoValues1);
        SimFluoValues2 = array2table(SimFluoValues2);
        SimFluoValues3 = array2table(SimFluoValues3);
%         
%         %add a column to specify the repression coefficient
        SimFluoValues1.Dose = zeros(height(SimFluoValues1),1);
        SimFluoValues2.Dose = zeros(height(SimFluoValues2),1);
        SimFluoValues3.Dose = zeros(height(SimFluoValues3),1);
%         
        SimFluoValues1.Repression_Coefficient = zeros(height(SimFluoValues1),1);
        SimFluoValues2.Repression_Coefficient = zeros(height(SimFluoValues2),1);
        SimFluoValues3.Repression_Coefficient = zeros(height(SimFluoValues3),1);
        
        SimFluoValues1.draw = zeros(height(SimFluoValues1),1);
        SimFluoValues2.draw = zeros(height(SimFluoValues2),1);
        SimFluoValues3.draw = zeros(height(SimFluoValues3),1);
% 
%         %set column names
        SimFluoValues1.Properties.VariableNames = {'SimFluoValues', 'Dose', 'Repression_Coefficient', 'Draw'};
        SimFluoValues2.Properties.VariableNames = {'SimFluoValues', 'Dose', 'Repression_Coefficient', 'Draw'};
        SimFluoValues3.Properties.VariableNames = {'SimFluoValues', 'Dose', 'Repression_Coefficient', 'Draw'};
        
        SimFluoValues1.dose = data1.dose;
        SimFluoValues2.dose = data2.dose;
        SimFluoValues3.dose = data3.dose;

        SimFluoValues1.Draw(:) = draw;
        SimFluoValues2.Draw(:) = draw;
        SimFluoValues3.Draw(:) = draw;
% 
%         %convert second column to string 
        SimFluoValues1 = convertvars(SimFluoValues1, "Repression_Coefficient", 'string');
        SimFluoValues2 = convertvars(SimFluoValues2, "Repression_Coefficient", 'string');
        SimFluoValues3 = convertvars(SimFluoValues3, "Repression_Coefficient", 'string');
        
        %insert the appropriate repression coefficient
        SimFluoValues1.Repression_Coefficient(:) = 'LacIrep1';
        SimFluoValues2.Repression_Coefficient(:) = 'LacIrep2';
        SimFluoValues3.Repression_Coefficient(:) = 'LacIrep3';
%         
%         %combine all tables to one
        if draw == 1
            SimFluoValues1combined = SimFluoValues1;
            SimFluoValues2combined = SimFluoValues2;
            SimFluoValues3combined = SimFluoValues3;
            SimFluoValuescombined = [SimFluoValues1; SimFluoValues2; SimFluoValues3];
        else
            SimFluoValues1combined = [SimFluoValues1combined; SimFluoValues1];
            SimFluoValues2combined = [SimFluoValues2combined; SimFluoValues2];
            SimFluoValues3combined = [SimFluoValues3combined; SimFluoValues3];
            SimFluoValuescombined = [SimFluoValuescombined; SimFluoValues1; SimFluoValues2; SimFluoValues3];
        end
% 
%         
%         
% 
%         %set the figure to be sceensize
%             set(figure(draw+3), 'Position', get(0, 'Screensize'));
%             %save figure
%             saveas(figure(draw+3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),num2str(draw),'Figure_draw_allrep', '.jpg']);

            counter = counter + 1
            fprintf('End run', counter)

end
%             %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Pesto/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep1_ec50_Pesto_statistics', '.mat'], 'SimFluoValues1combined');
%             writetable(SimFluoValues1combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep1_ec50', '.csv']);
%              %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Pesto/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep2_ec50_Pesto_statistics', '.mat'], 'SimFluoValues2combined');
%             writetable(SimFluoValues2combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep2_ec50', '.csv']);
%             %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Pesto/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep3_ec50_Pesto_statistics', '.mat'], 'SimFluoValues3combined');
%             writetable(SimFluoValues3combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep3_ec50', '.csv']);
%            
%             save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'SimulationFluoValues_allrep_ec50', '.mat'], 'SimFluoValuescombined');
           

%% Calculate ec50 for pesto data

counter = 0;
for k = 1:num_draws
    x = SimFluoValues1combined.dose(SimFluoValues1combined.Draw == k);
    y = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_1_pesto(k,1) = xq(:,closest_index);
    ec50_1_pesto(k,2) = log10(ec50_1_pesto(k,1));
    ec50_1_pesto(k,3) = mid_fluo;
    ec50_1_pesto(k,4) = interpol(closest_index);
    ec50_1_pesto(k,5) = (abs(ec50_1_pesto(k,3) - ec50_1_pesto(k,4))) / mid_fluo;
    ec50_1_pesto(k,6) = k;
    counter = counter +1
end
counter = 0;
for k = 1:num_draws
    x = SimFluoValues2combined.dose(SimFluoValues2combined.Draw == k);
    y = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_2_pesto(k,1) = xq(:,closest_index);
    ec50_2_pesto(k,2) = log10(ec50_2_pesto(k,1));
    ec50_2_pesto(k,3) = mid_fluo;
    ec50_2_pesto(k,4) = interpol(closest_index);
    ec50_2_pesto(k,5) = (abs(ec50_2_pesto(k,3) - ec50_2_pesto(k,4))) / mid_fluo;
    ec50_2_pesto(k,6) = k;
    counter = counter +1
end
counter = 0;

for k = 1:num_draws
    x = SimFluoValues3combined.dose(SimFluoValues3combined.Draw == k);
    y = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 5 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_3_pesto(k,1) = xq(:,closest_index);
    ec50_3_pesto(k,2) = log10(ec50_3_pesto(k,1));
    ec50_3_pesto(k,3) = mid_fluo;
    ec50_3_pesto(k,4) = interpol(closest_index);
    ec50_3_pesto(k,5) = (abs(ec50_3_pesto(k,3) - ec50_3_pesto(k,4))) / mid_fluo;
    ec50_3_pesto(k,6) = k;
    counter = counter +1
end

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Pesto/', datestr(now, 'dd-mmm-yyyy'),'ec50_1_random_pesto_statistics', '.mat'], 'ec50_1_pesto');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Pesto/', datestr(now, 'dd-mmm-yyyy'),'ec50_2_random_pesto_statistics', '.mat'], 'ec50_2_pesto');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Pesto/', datestr(now, 'dd-mmm-yyyy'),'ec50_3_random_pesto_statistics', '.mat'], 'ec50_3_pesto');


%% Statistics

%standard error
SEM50_1_pesto = std(ec50_1_pesto(:,2))/sqrt(length(ec50_1_pesto(:,2)));
SEM50_2_pesto = std(ec50_2_pesto(:,2))/sqrt(length(ec50_2_pesto(:,2)));
SEM50_3_pesto = std(ec50_3_pesto(:,2))/sqrt(length(ec50_3_pesto(:,2)));


%T-Score
ts50_1_pesto = tinv([0.05 0.95], length(ec50_1_pesto(:,2))-1);
ts50_2_pesto = tinv([0.05 0.95], length(ec50_2_pesto(:,2))-1);
ts50_3_pesto = tinv([0.05 0.95], length(ec50_3_pesto(:,2))-1);

%Confidence Interval
CI50_1_pesto = mean(ec50_1_pesto(:,2)) + ts50_1_pesto*SEM50_1_pesto;
CI50_2_pesto = mean(ec50_2_pesto(:,2)) + ts50_2_pesto*SEM50_2_pesto;
CI50_3_pesto = mean(ec50_3_pesto(:,2)) + ts50_3_pesto*SEM50_3_pesto;

%Standard deviation
std50_1_pesto = std(ec50_1_pesto(:,2));
std50_2_pesto = std(ec50_2_pesto(:,2));
std50_3_pesto = std(ec50_3_pesto(:,2));

%Interquantile range
r50_1_pesto = iqr(ec50_1_pesto(:,2));
r50_2_pesto = iqr(ec50_2_pesto(:,2));
r50_3_pesto = iqr(ec50_3_pesto(:,2));

%Mean
mean_ec50_1_pesto = mean(ec50_1_pesto(:,2));
mean_ec50_2_pesto = mean(ec50_2_pesto(:,2));
mean_ec50_3_pesto = mean(ec50_3_pesto(:,2));

%Variance
var_ec50_1_pesto = var(ec50_1_pesto(:,2));
var_ec50_2_pesto = var(ec50_2_pesto(:,2));
var_ec50_3_pesto = var(ec50_3_pesto(:,2));

%ANOVA
all_ec50s_pesto = [ec50_1_pesto(:,2) ec50_2_pesto(:,2) ec50_3_pesto(:,2)];
group = ["rep1" "rep2" "rep3"];
[p_50_pesto,tbl_50_pesto,stats_50_pesto] = anova1(all_ec50s_pesto, group);

ec50_1_2_pesto = [ec50_1_pesto(:,2) ec50_2_pesto(:,2)];
group = ["rep1" "rep2"];
[p_50_1_2_pesto, tbl_50_1_2_pesto, stats_50_1_2_pesto] = anova1(ec50_1_2_pesto, group);

ec50_1_3_pesto = [ec50_1_pesto(:,2) ec50_3_pesto(:,2)];
group = ["rep1" "rep3"];
[p_50_1_3_pesto, tbl_50_1_3_pesto, stats_50_1_3_pesto] = anova1(ec50_1_3_pesto, group);

ec50_2_3_pesto = [ec50_2_pesto(:,2) ec50_3_pesto(:,2)];
group = ["rep2" "rep3"];
[p_50_2_3_pesto, tbl_50_2_3_pesto, stats_50_2_3_pesto] = anova1(ec50_2_3_pesto, group);

%Coefficient of Variation
CV50_1_pesto = (std(ec50_1_pesto(:,2))/mean(ec50_1_pesto(:,2))) * 100;
CV50_2_pesto = (std(ec50_2_pesto(:,2))/mean(ec50_2_pesto(:,2))) * 100;
CV50_3_pesto = (std(ec50_3_pesto(:,2))/mean(ec50_3_pesto(:,2))) * 100;

%Histograms
subplot(1,3,1)
hist_ec50_1_pesto = histogram(ec50_1_pesto(:,2));
title('EC50 repression coefficient 1 - randomly drawn parameters from PESTO data')
subplot(1,3,2)
hist_ec50_2_pesto = histogram(ec50_2_pesto(:,2));
title('EC50 repression coefficient 2 - randomly drawn parameters from PESTO data')
subplot(1,3,3)
hist_ec50_3_pesto = histogram(ec50_3_pesto(:,2));
title('EC50 repression coefficient 3 - randomly drawn parameters from PESTO data')

%save in one structure

statistics_PESTO.sem1 = SEM50_1_pesto;
statistics_PESTO.sem2 = SEM50_2_pesto;
statistics_PESTO.sem3 = SEM50_3_pesto;

statistics_PESTO.tscore1 = ts50_1_pesto;
statistics_PESTO.tscore2 = ts50_2_pesto;
statistics_PESTO.tscore3 = ts50_3_pesto;

statistics_PESTO.CI1 = CI50_1_pesto;
statistics_PESTO.CI2 = CI50_2_pesto;
statistics_PESTO.CI3 = CI50_3_pesto;

statistics_PESTO.std1 = std50_1_pesto;
statistics_PESTO.std2 = std50_2_pesto;
statistics_PESTO.std3 = std50_3_pesto;

statistics_PESTO.interquartile_range1 = r50_1_pesto;
statistics_PESTO.interquartile_range2 = r50_2_pesto;
statistics_PESTO.interquartile_range3 = r50_3_pesto;

statistics_PESTO.mean1 = mean_ec50_1_pesto;
statistics_PESTO.mean2 = mean_ec50_2_pesto;
statistics_PESTO.mean3 = mean_ec50_3_pesto;

statistics_PESTO.var1 = var_ec50_1_pesto;
statistics_PESTO.var2 = var_ec50_2_pesto;
statistics_PESTO.var3 = var_ec50_3_pesto;

statistics_PESTO.ANOVA.all.p = p_50_pesto;
statistics_PESTO.ANOVA.all.tbl = tbl_50_pesto;
statistics_PESTO.ANOVA.all.stats = stats_50_pesto;

statistics_PESTO.ANOVA.rep1_rep3.p = p_50_1_3_pesto;
statistics_PESTO.ANOVA.rep1_rep3.tbl = tbl_50_1_3_pesto;
statistics_PESTO.ANOVA.rep1_rep3.stats = stats_50_1_3_pesto;

statistics_PESTO.ANOVA.rep1_rep2.p = p_50_1_2_pesto;
statistics_PESTO.ANOVA.rep1_rep2.tbl = tbl_50_1_2_pesto;
statistics_PESTO.ANOVA.rep1_rep2.stats = stats_50_1_2_pesto;

statistics_PESTO.ANOVA.rep2_rep3.p = p_50_2_3_pesto;
statistics_PESTO.ANOVA.rep2_rep3.tbl = tbl_50_2_3_pesto;
statistics_PESTO.ANOVA.rep2_rep3.stats = stats_50_2_3_pesto;

statistics_PESTO.CV1 = CV50_1_pesto;
statistics_PESTO.CV2 = CV50_2_pesto;
statistics_PESTO.CV3 = CV50_3_pesto;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Statistics/Parameters from PESTO/',datestr(now, 'dd-mmm-yyyy'),'_statistics_PESTO.mat'],'statistics_PESTO')


