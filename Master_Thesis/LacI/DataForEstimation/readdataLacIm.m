function data_LacI = readdataLacI

% used to generate the one .mat file per experiment (strains measured on
% the same day) from the excel or facs data

data.time = 20; % timepoint(s) of measurement

data_raw = readtable('/DataForEstimation/Citrine_size_corrected.csv');

data.tdh3 = table2array(data_raw(86,"Mean"));
data.empty = table2array(data_raw(end, "Mean"));
data.means = data_raw(1:12,3);
data.std = data_raw(1:12,6);

modulenames{1} = 'PA_LacI-P4Lacn_cit';

dose = data_raw(1:12,9);
timepoint = data_raw(1:12,13);

save('/Users/claudiobroger/Documents/ETH/Master_Thesis/LacI/DataForEstimation/data.mat', 'data', 'modulenames', 'dose', 'timepoint')
end





