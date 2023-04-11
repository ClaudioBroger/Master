path = '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/';
dataPos = strcat(path, "13-Mar-2023SimulationFluoValues_rep1_.mat");
dataPos2 = strcat(path, "13-Mar-2023SimulationFluoValues_rep2_.mat");
dataPos3 = strcat(path, "13-Mar-2023SimulationFluoValues_rep3_.mat");
load(dataPos);
load(dataPos2);
load(dataPos3);
num_draws = 20;
counter = 0;


%% EC50 for randomly drawn parameter sets
counter = 0;
for k = 1:num_draws
    x = SimFluoValues1combined.dose(SimFluoValues1combined.Draw == k);
    y = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_1(k,1) = xq(:,closest_index);
    ec50_1(k,2) = log10(ec50_1(k,1));
    ec50_1(k,3) = mid_fluo;
    ec50_1(k,4) = interpol(closest_index);
    ec50_1(k,5) = (abs(ec50_1(k,3) - ec50_1(k,4))) / mid_fluo;
    ec50_1(k,6) = k;
    counter = counter +1
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/', datestr(now, 'dd-mmm-yyyy'),'ec50_1_random_200', '.mat'], 'ec50_1');

counter = 0;
for k = 1:num_draws
    x = SimFluoValues2combined.dose(SimFluoValues2combined.Draw == k);
    y = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_2(k,1) = xq(:,closest_index);
    ec50_2(k,2) = log10(ec50_2(k,1));
    ec50_2(k,3) = mid_fluo;
    ec50_2(k,4) = interpol(closest_index);
    ec50_2(k,5) = (abs(ec50_2(k,3) - ec50_2(k,4))) / mid_fluo;
    ec50_2(k,6) = k;
    counter = counter +1
end

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/', datestr(now, 'dd-mmm-yyyy'),'ec50_2_random_200', '.mat'], 'ec50_2');

counter = 0;
for k = 1:num_draws
    x = SimFluoValues3combined.dose(SimFluoValues3combined.Draw == k);
    y = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 5 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_3(k,1) = xq(:,closest_index);
    ec50_3(k,2) = log10(ec50_3(k,1));
    ec50_3(k,3) = mid_fluo;
    ec50_3(k,4) = interpol(closest_index);
    ec50_3(k,5) = (abs(ec50_3(k,3) - ec50_3(k,4))) / mid_fluo;
    ec50_3(k,6) = k;
    counter = counter +1
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/', datestr(now, 'dd-mmm-yyyy'),'ec50_3_random_200', '.mat'], 'ec50_3');


C = linspecer(num_draws);
figure(1)
hold on
for k = 1:num_draws
    plot(ec50_1(k,2), ec50_1(k,4), '-k*','Color', C(k,:),'DisplayName',strcat('EC50-', num2str(k)),'MarkerSize',12)
end


figure(2)
hold on
for k = 1:num_draws
    plot(ec50_2(k,2), ec50_2(k,4), '-k*','Color', C(k,:),'DisplayName',strcat('EC50-', num2str(k)),'MarkerSize',12)
end

figure(3)
hold on
for k = 1:num_draws
    plot(ec50_3(k,2), ec50_3(k,4), '-k*','Color', C(k,:),'DisplayName',strcat('EC50-', num2str(k)),'MarkerSize',12)
end

%% EC90 for randomly drawn parameter sets

for k = 1:num_draws
    x = SimFluoValues1combined.dose(SimFluoValues1combined.Draw == k);
    y = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.9);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec90_1(k,1) = xq(:,closest_index);
    ec90_1(k,2) = log10(ec90_1(k,1));
    ec90_1(k,3) = mid_fluo;
    ec90_1(k,4) = interpol(closest_index);
    ec90_1(k,5) = (abs(ec90_1(k,3) - ec90_1(k,4))) / mid_fluo;
    ec90_1(k,6) = k;
    counter = counter +1
end

for k = 1:num_draws
    x = SimFluoValues2combined.dose(SimFluoValues2combined.Draw == k);
    y = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.9);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec90_2(k,1) = xq(:,closest_index);
    ec90_2(k,2) = log10(ec90_2(k,1));
    ec90_2(k,3) = mid_fluo;
    ec90_2(k,4) = interpol(closest_index);
    ec90_2(k,5) = (abs(ec90_2(k,3) - ec90_2(k,4))) / mid_fluo;
    ec90_2(k,6) = k;
    counter = counter +1
end


for k = 1:num_draws
    x = SimFluoValues3combined.dose(SimFluoValues3combined.Draw == k);
    y = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.9);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec90_3(k,1) = xq(:,closest_index);
    ec90_3(k,2) = log10(ec90_3(k,1));
    ec90_3(k,3) = mid_fluo;
    ec90_3(k,4) = interpol(closest_index);
    ec90_3(k,5) = (abs(ec90_3(k,3) - ec90_3(k,4))) / mid_fluo;
    ec90_3(k,6) = k;
    counter = counter +1
end

%% EC50 for parameter sets drawn from 10% below and 10% above optimal parameter set

for k = 1:num_draws
    x = SimFluoValues1combined.Dose(SimFluoValues1combined.Draw == k);
    y = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_1_around_paraopt(k,1) = xq(:,closest_index);
    ec50_1_around_paraopt(k,2) = log10(ec50_1_around_paraopt(k,1));
    ec50_1_around_paraopt(k,3) = mid_fluo;
    ec50_1_around_paraopt(k,4) = interpol(closest_index);
    ec50_1_around_paraopt(k,5) = (abs(ec50_1_around_paraopt(k,3) - ec50_1_around_paraopt(k,4))) / mid_fluo;
end

for k = 1:num_draws
    x = SimFluoValues2combined.Dose(SimFluoValues2combined.Draw == k);
    y = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_2_around_paraopt(k,1) = xq(:,closest_index);
    ec50_2_around_paraopt(k,2) = log10(ec50_2_around_paraopt(k,1));
    ec50_2_around_paraopt(k,3) = mid_fluo;
    ec50_2_around_paraopt(k,4) = interpol(closest_index);
    ec50_2_around_paraopt(k,5) = (abs(ec50_2_around_paraopt(k,3) - ec50_2_around_paraopt(k,4))) / mid_fluo;
end


for k = 1:num_draws
    x = SimFluoValues3combined.Dose(SimFluoValues3combined.Draw == k);
    y = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_3_around_paraopt(k,1) = xq(:,closest_index);
    ec50_3_around_paraopt(k,2) = log10(ec50_3_around_paraopt(k,1));
    ec50_3_around_paraopt(k,3) = mid_fluo;
    ec50_3_around_paraopt(k,4) = interpol(closest_index);
    ec50_3_around_paraopt(k,5) = (abs(ec50_3_around_paraopt(k,3) - ec50_3_around_paraopt(k,4))) / mid_fluo;
end


x = SimFluoValues1_paraopt.dose;
y = SimFluoValues1_paraopt.SimFluoValues;
mid_fluo = abs(min(y) + max(y) *0.5);
xq = (min(x) : 1 : max(x));
interpol = interp1(x,y,xq);
[minvalue closest_index] = min(abs(interpol - mid_fluo));
ec50_1_paraopt(1,1) = xq(:,closest_index);
ec50_1_paraopt(1,2) = log10(ec50_1_paraopt(1,1));
ec50_1_paraopt(1,3) = mid_fluo;
ec50_1_paraopt(1,4) = interpol(closest_index);
ec50_1_paraopt(1,5) = (abs(ec50_1_paraopt(1,3) - ec50_1_paraopt(1,4))) / mid_fluo;


x = SimFluoValues2_paraopt.dose;
y = SimFluoValues2_paraopt.SimFluoValues;
mid_fluo = abs(min(y) + max(y) *0.5);
xq = (min(x) : 1 : max(x));
interpol = interp1(x,y,xq);
[minvalue closest_index] = min(abs(interpol - mid_fluo));
ec50_2_paraopt(1,1) = xq(:,closest_index);
ec50_2_paraopt(1,2) = log10(ec50_2_paraopt(1,1));
ec50_2_paraopt(1,3) = mid_fluo;
ec50_2_paraopt(1,4) = interpol(closest_index);
ec50_2_paraopt(1,5) = (abs(ec50_2_paraopt(1,3) - ec50_2_paraopt(1,4))) / mid_fluo;


x = SimFluoValues3_paraopt.dose;
y = SimFluoValues3_paraopt.SimFluoValues;
mid_fluo = abs(min(y) + max(y) *0.5);
xq = (min(x) : 1 : max(x));
interpol = interp1(x,y,xq);
[minvalue closest_index] = min(abs(interpol - mid_fluo));
ec50_3_paraopt(1,1) = xq(:,closest_index);
ec50_3_paraopt(1,2) = log10(ec50_3_paraopt(1,1));
ec50_3_paraopt(1,3) = mid_fluo;
ec50_3_paraopt(1,4) = interpol(closest_index);
ec50_3_paraopt(1,5) = (abs(ec50_3_paraopt(1,3) - ec50_3_paraopt(1,4))) / mid_fluo;

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/28. March/', datestr(now, 'dd-mmm-yyyy'),'ec50_1_random', '.mat'], 'ec50_1');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/28. March/', datestr(now, 'dd-mmm-yyyy'),'ec50_2_random', '.mat'], 'ec50_2');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/28. March/', datestr(now, 'dd-mmm-yyyy'),'ec50_3_random', '.mat'], 'ec50_3');

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/22. March/', datestr(now, 'dd-mmm-yyyy'),'ec90_1_random', '.mat'], 'ec90_1');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/22. March/', datestr(now, 'dd-mmm-yyyy'),'ec90_2_random', '.mat'], 'ec90_2');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/22. March/', datestr(now, 'dd-mmm-yyyy'),'ec90_3_random', '.mat'], 'ec90_3');
   


