
path = '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/';
dataPos = strcat(path, "13-Mar-2023SimulationFluoValues_rep1_.mat");
dataPos2 = strcat(path, "13-Mar-2023SimulationFluoValues_rep2_.mat");
dataPos3 = strcat(path, "13-Mar-2023SimulationFluoValues_rep3_.mat");
load(dataPos);
load(dataPos2);
load(dataPos3);
num_draws = 20;

%Comparing maximum deviations between response curves with parameters
%randomly drawn


for k = 1:num_draws
    for j = 1:num_draws
        maximum_deviations1(j,k) = max(abs(SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k) - SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == j)));
        maximum_deviations2(j,k) = max(abs(SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k) - SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == j)));
        maximum_deviations3(j,k) = max(abs(SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k) - SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == j)));
    end
end
m1 = tril(maximum_deviations1);
m2 = tril(maximum_deviations2);
m3 = tril(maximum_deviations3);

for k = 1:num_draws
    deviates1(k) = sum(m1(k,:)) + sum(m1(:,k));
    deviates2(k) = sum(m2(k,:)) + sum(m2(:,k));
    deviates3(k) = sum(m3(k,:)) + sum(m3(:,k));
end
[most_deviate1, Index_most_deviate1] = max(deviates1);
[min_deviate1, Index_min_deviate1] = min(deviates1);
[most_deviate2, Index_most_deviate2] = max(deviates2);
[min_deviate2, Index_min_deviate2] = min(deviates2);
[most_deviate3, Index_most_deviate3] = max(deviates3);
[min_deviate3, Index_min_deviate3] = min(deviates3);




%Comparing maximum deviations between response curves with parameters
%randomly drawn and the experimental data

dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');
dataPos = strcat(dataPath, "data.mat");
data = load(dataPos);

for k = 1:num_draws
    for j = 1:num_draws
        maximum_deviations_from_data1(j,k) = max(abs(SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k) - data.data.means));
        maximum_deviations_from_data2(j,k) = max(abs(SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k) - data.data.means));
        maximum_deviations_from_data3(j,k) = max(abs(SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k) - data.data.means));
    end
end
m_from_data1 = tril(maximum_deviations_from_data1);
m_from_data2 = tril(maximum_deviations_from_data2);
m_from_data3 = tril(maximum_deviations_from_data3);

for k = 1:num_draws
    deviates_from_data1(k) = sum(m_from_data1(k,:)) + sum(m_from_data1(:,k));
    deviates_from_data2(k) = sum(m_from_data2(k,:)) + sum(m_from_data2(:,k));
    deviates_from_data3(k) = sum(m_from_data3(k,:)) + sum(m_from_data3(:,k));
end
[most_deviate_from_data1, Index_most_deviate_from_data1] = max(deviates_from_data1);
[min_deviate_from_data1, Index_min_deviate_from_data1] = min(deviates_from_data1);
[most_deviate_from_data2, Index_most_deviate_from_data2] = max(deviates_from_data2);
[min_deviate_from_data2, Index_min_deviate_from_data2] = min(deviates_from_data2);
[most_deviate_from_data3, Index_most_deviate_from_data3] = max(deviates_from_data3);
[min_deviate_from_data3, Index_min_deviate_from_data3] = min(deviates_from_data3);


path = '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/';
dataPos = strcat(path, "13-Mar-2023_rand_parameter_uncertainty_analysis.mat");
load(dataPos)

rand_parameter = table2array(rand_parameter);
for k = 1:num_draws
    for j = 1:width(rand_parameter)
        differences_parameter_to_min_deviation1(k,j) = abs(rand_parameter(Index_min_deviate1,j) - rand_parameter(k,j));
    end
end
