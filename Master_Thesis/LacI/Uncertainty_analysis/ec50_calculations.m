path = '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/';
dataPos = strcat(path, "13-Mar-2023SimulationFluoValues_rep1_.mat");
dataPos2 = strcat(path, "13-Mar-2023SimulationFluoValues_rep2_.mat");
dataPos3 = strcat(path, "13-Mar-2023SimulationFluoValues_rep3_.mat");
load(dataPos);
load(dataPos2);
load(dataPos3);
num_draws = 20;



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
end

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
end


for k = 1:num_draws
    x = SimFluoValues3combined.dose(SimFluoValues3combined.Draw == k);
    y = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_3(k,1) = xq(:,closest_index);
    ec50_3(k,2) = log10(ec50_3(k,1));
    ec50_3(k,3) = mid_fluo;
    ec50_3(k,4) = interpol(closest_index);
    ec50_3(k,5) = (abs(ec50_3(k,3) - ec50_3(k,4))) / mid_fluo;
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




