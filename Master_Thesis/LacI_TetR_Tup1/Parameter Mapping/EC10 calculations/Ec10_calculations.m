path = '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC50 calculations';
dataPos = strcat(path, "27-Apr-2023SimFluoValues1_LacI_TetRTup1_EC50.mat");
dataPos2 = strcat(path, "27-Apr-2023SimFluoValues2_LacI_TetRTup1_EC50.mat");
dataPos3 = strcat(path, "27-Apr-2023SimFluoValues3_LacI_TetRTup1_EC50.mat");
load(dataPos);
load(dataPos2);
load(dataPos3);
num_draws = 200;
counter = 0;
%% EC50 for randomly drawn parameter sets
counter = 0;
for k = 1:num_draws
    x = SimFluoValues1combined.Dose(SimFluoValues1combined.Draw == k);
    y = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.1);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec10_1(k,1) = xq(:,closest_index);
    ec10_1(k,2) = log10(ec10_1(k,1));
    ec10_1(k,3) = mid_fluo;
    ec10_1(k,4) = interpol(closest_index);
    ec10_1(k,5) = (abs(ec10_1(k,3) - ec10_1(k,4))) / mid_fluo;
    ec10_1(k,6) = k;
    counter = counter +1
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ec10_1_TetRTup1', '.mat'], 'ec10_1');

counter = 0;
for k = 1:num_draws
    x = SimFluoValues2combined.Dose(SimFluoValues2combined.Draw == k);
    y = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec10_2(k,1) = xq(:,closest_index);
    ec10_2(k,2) = log10(ec10_2(k,1));
    ec10_2(k,3) = mid_fluo;
    ec10_2(k,4) = interpol(closest_index);
    ec10_2(k,5) = (abs(ec10_2(k,3) - ec10_2(k,4))) / mid_fluo;
    ec10_2(k,6) = k;
    counter = counter +1
end

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ec10_2_TetRTup1', '.mat'], 'ec10_2');

counter = 0;
for k = 1:num_draws
    x = SimFluoValues3combined.Dose(SimFluoValues3combined.Draw == k);
    y = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec10_3(k,1) = xq(:,closest_index);
    ec10_3(k,2) = log10(ec10_3(k,1));
    ec10_3(k,3) = mid_fluo;
    ec10_3(k,4) = interpol(closest_index);
    ec10_3(k,5) = (abs(ec10_3(k,3) - ec10_3(k,4))) / mid_fluo;
    ec10_3(k,6) = k;
    counter = counter +1
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/EC10 calculations/', datestr(now, 'dd-mmm-yyyy'),'ec10_3_TetRTup1', '.mat'], 'ec10_3');
