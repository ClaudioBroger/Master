path = '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/';
dataPos = strcat(path, "13-Mar-2023SimulationFluoValues_rep1_.mat");
dataPos2 = strcat(path, "13-Mar-2023SimulationFluoValues_rep2_.mat");
dataPos3 = strcat(path, "13-Mar-2023SimulationFluoValues_rep3_.mat");
load(dataPos);
load(dataPos2);
load(dataPos3);
num_draws = 20;

for k= 1:num_draws
        [hill(k,1) ec50(k,1)] = doseResponse(SimFluoValues1combined.Dose(SimFluoValues1combined.Draw == k), SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k));

end

for k= 1:num_draws
        [hill2(k,1) ec502(k,1)] = doseResponse(SimFluoValues2combined.Dose(SimFluoValues2combined.Draw == k), SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k));

end

for k= 1:num_draws
        [hill3(k,1) ec503(k,1)] = doseResponse(SimFluoValues3combined.Dose(SimFluoValues3combined.Draw == k), SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k));

end


plot(ec50)
plot(hill)
plot(ec502)
plot(hill2)
plot(ec503)
plot(hill3)




