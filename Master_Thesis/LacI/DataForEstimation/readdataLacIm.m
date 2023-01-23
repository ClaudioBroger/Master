function data_LacI = readdataLacI

% used to generate the one .mat file per experiment (strains measured on
% the same day) from the excel or facs data

data.time = 20; % timepoint(s) of measurement

data_raw = readtable('/DataForEstimation/Citrine_size_corrected.csv');

data.tdh3 = table2array(data_raw(86,"Mean"));
data.empty = table2array(data_raw(end, "Mean"));
data.means = table2array(data_raw(1:12,3));
data.std = table2array(data_raw(1:12,6));

modulenames{1} = 'PA_LacI-P4Lacn_cit';

dose = table2array(data_raw(1:12,9))*1000;
timepoint = table2array(data_raw(1:12,13));

save('LacI/DataForEstimation/data.mat', 'data', 'modulenames', 'dose', 'timepoint')

data_W220F.time = 20;

data_W220F.tdh3 = table2array(data_raw(86,"Mean"));
data_W220F.empty = table2array(data_raw(end, "Mean"));
data_W220F.means = table2array(data_raw(13:24,3));
data_W220F.std = table2array(data_raw(13:24,6));

modulenames_W220F{1} = 'PAct1_LacI(W220F)_tCyc1';

dose_W220F = table2array(data_raw(13:24,9))*1000;
timepoint_W220F = table2array(data_raw(13:24,13));

save('LacI/DataForEstimation/data_W220F.mat', 'data_W220F', 'modulenames_W220F', 'dose_W220F', 'timepoint_W220F')

data_W220F_Q60G_T167A.time = 20;

data_W220F_Q60G_T167A.tdh3 = table2array(data_raw(86,"Mean"));
data_W220F_Q60G_T167A.empty = table2array(data_raw(end, "Mean"));
data_W220F_Q60G_T167A.means = table2array(data_raw(25:36,3));
data_W220F_Q60G_T167A.std = table2array(data_raw(25:36,6));

modulenames_W220F_Q60G_T167A{1} = 'P4Lacn.2_LacI(W220F,Q60G, T167A)_tCyc1'

dose_W220F_Q60G_T167A = table2array(data_raw(25:36,9))*1000;
timepoint_W220F_Q60G_T167A = table2array(data_raw(25:36,13));

save('LacI/DataForEstimation/data_W220F_Q60G_T167A.mat', 'data_W220F_Q60G_T167A', 'modulenames_W220F_Q60G_T167A', 'dose_W220F_Q60G_T167A', 'timepoint_W220F_Q60G_T167A')

data_pt7.time = 20;

data_pt7.tdh3 = table2array(data_raw(86,"Mean"));
data_pt7.empty = table2array(data_raw(end, "Mean"));
data_pt7.means = table2array(data_raw(37:48,3));
data_pt7.std = table2array(data_raw(37:48,6));

modulenames_pt7{1} = 'PAct1_citrine-pt7-LacI(W220F, Q60G, T167A)'

dose_pt7 = table2array(data_raw(37:48,9))*1000;
timepoint_pt7 = table2array(data_raw(37:48,13));

save('LacI/DataForEstimation/data_pt7.mat', 'data_pt7', 'modulenames_pt7', 'dose_pt7', 'timepoint_pt7')

end



