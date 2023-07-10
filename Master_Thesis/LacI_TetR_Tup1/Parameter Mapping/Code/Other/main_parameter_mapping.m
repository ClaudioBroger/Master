clear all; clc;

%Reset, reload and compile model
sbioreset;
simple_model_LacI_TetR;
sbiosaveproject('simple_model');
model = sbioloadproject('simple_model');

sbioaccelerate(model.mw_sbmod1)

%Modelsettings
ModelSettings_LacI_TetR;

%draw parameter values
draw_parameters_LacI_TetR_parameter_mapping;

%Simulate with drawn parameters
Simulate_random_drawn_parameters_parameter_mapping;


counter = 0;
for k = 1:num_draws
    x = SimFluoValues1combined.Dose(SimFluoValues1combined.Draw == k);
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
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/', datestr(now, 'dd-mmm-yyyy'),'ec50_1_30', '.mat'], 'ec50_1');

counter = 0;
for k = 1:num_draws
    x = SimFluoValues2combined.Dose(SimFluoValues2combined.Draw == k);
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

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/', datestr(now, 'dd-mmm-yyyy'),'ec50_2_30', '.mat'], 'ec50_2');

counter = 0;
for k = 1:num_draws
    x = SimFluoValues3combined.Dose(SimFluoValues3combined.Draw == k);
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
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/', datestr(now, 'dd-mmm-yyyy'),'ec50_3_30', '.mat'], 'ec50_3');


C = linspecer(num_draws);
figure(2)
hold on
for k = 1:num_draws
    plot(ec50_1(k,2), ec50_1(k,4), '-k*','Color', C(k,:),'DisplayName',strcat('EC50-', num2str(k)),'MarkerSize',12)
end


figure(3)
hold on
for k = 1:num_draws
    plot(ec50_2(k,2), ec50_2(k,4), '-k*','Color', C(k,:),'DisplayName',strcat('EC50-', num2str(k)),'MarkerSize',12)
end

figure(4)
hold on
for k = 1:num_draws
    plot(ec50_3(k,2), ec50_3(k,4), '-k*','Color', C(k,:),'DisplayName',strcat('EC50-', num2str(k)),'MarkerSize',12)
end

