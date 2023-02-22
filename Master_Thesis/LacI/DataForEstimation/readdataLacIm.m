function data_LacI = readdataLacI

% used to generate the one .mat file per experiment (strains measured on
% the same day) from the excel or facs data

data.time = 1200; % timepoint(s) of measurement

data_raw = readtable('/DataForEstimation/Citrine_size_corrected.csv');

data.tdh3 = table2array(data_raw(86,"Mean"));
data.empty = table2array(data_raw(end, "Mean"));
data.means = data_raw.Mean(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_LacI'));
data.std = data_raw.Stdev(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_LacI'));

modulenames{1} = 'PA_LacI-P4Lacn_cit';

dose = data_raw.Dose(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_LacI'))*1000;
timepoint = data_raw.Time(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_LacI'));

save('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/data.mat','data', 'modulenames', 'dose', 'timepoint')

%data_W220F.tdh3 = data.tdh3;
%data_W220F.empty = data.tdh3;
data.means = data_raw.Mean(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_LacI(W220F)'));
data.std = data_raw.Stdev(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_LacI(W220F)'));

modulenames{1} = 'PAct1_LacI(W220F)_tCyc1';

dose = data_raw.Dose(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_LacI(W220F)'))*1000;
timepoint = data_raw.Time(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_LacI(W220F)'));

save('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/data_W220F.mat','data', 'modulenames', 'dose', 'timepoint')


%data_W220F_Q60G_T167A.time = data.time;

%data_W220F_Q60G_T167A.tdh3 = data.tdh3;
%data_W220F_Q60G_T167A.empty = data.empty;
data.means = data_raw.Mean(strcmp(data_raw.Strain, 'P4Lacn.2_cit + P4Lacn.2_LacI(W220F,Q60G, T167A)'));
data.std = data_raw.Stdev(strcmp(data_raw.Strain, 'P4Lacn.2_cit + P4Lacn.2_LacI(W220F,Q60G, T167A)'));

modulenames{1} = 'P4Lacn.2_LacI(W220F,Q60G, T167A)_tCyc1';

dose = data_raw.Dose(strcmp(data_raw.Strain, 'P4Lacn.2_cit + P4Lacn.2_LacI(W220F,Q60G, T167A)'))*1000;
timepoint = data_raw.Time(strcmp(data_raw.Strain, 'P4Lacn.2_cit + P4Lacn.2_LacI(W220F,Q60G, T167A)'));

save('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/data_W220F_Q60G_T167A.mat','data', 'modulenames', 'dose', 'timepoint')

%data_pt7.time = data.time;

%data_pt7.tdh3 = data.tdh3;
%data_pt7.empty = data.empty;
data.means = data_raw.Mean(strcmp(data_raw.Strain, 'P3Lacn.5_cit + PAct1_citrine-pt7-LacI(W220F, Q60G, T167A)'));
data.std = data_raw.Stdev(strcmp(data_raw.Strain, 'P3Lacn.5_cit + PAct1_citrine-pt7-LacI(W220F, Q60G, T167A)'));

modulenames{1} = 'PAct1_citrine-pt7-LacI(W220F, Q60G, T167A)';

dose = data_raw.Dose(strcmp(data_raw.Strain, 'P3Lacn.5_cit + PAct1_citrine-pt7-LacI(W220F, Q60G, T167A)'))*1000;
timepoint = data_raw.Time(strcmp(data_raw.Strain, 'P3Lacn.5_cit + PAct1_citrine-pt7-LacI(W220F, Q60G, T167A)'));

save('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/data_pt7.mat','data', 'modulenames', 'dose', 'timepoint')

%data_pt7_5circuit.time = data.time;

%data_pt7_5circuit.tdh3 = data.tdh3;
%data_pt7_5circuit.empty = data.empty;
data.means = data_raw.Mean(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_citrine-pt7-LacI(W220F, Q60G, T167A)'));
data.std = data_raw.Stdev(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_citrine-pt7-LacI(W220F, Q60G, T167A)'));

modulenames{1} = 'P4Lacn.2_citrine-pt7-LacI(W220F, Q60G, T167A)';

dose = data_raw.Dose(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_citrine-pt7-LacI(W220F, Q60G, T167A)'))*1000;
timepoint = data_raw.Time(strcmp(data_raw.Strain, 'P4Lacn.2_cit + PAct1_citrine-pt7-LacI(W220F, Q60G, T167A)'));

save('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/data_pt7_5circuit.mat','data', 'modulenames', 'dose', 'timepoint')

end



