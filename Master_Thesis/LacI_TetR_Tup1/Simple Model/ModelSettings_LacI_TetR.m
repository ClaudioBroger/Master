%%% Model & Estimation Settings %%%

Settings = [];

%General
Settings.model.modelproj = 'model.mw_sbmod1';
[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
[ParaValues] = cell2mat(get(model.mw_sbmod1.parameters, {'Value'}));
ParametersTable = readtable('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Simple Model/parameters_0_LacI_TetR.xlsx');

% Specific
[~,Locb] = ismember(ParaNames, ParametersTable.names); %tells you how paramspecs should be changed
paramSpecs = ParametersTable(Locb, :);
Settings.model.PIdx = find(strcmp(paramSpecs.estimate, 'yes'));
Settings.model.PIdxfixed = find(strcmp(paramSpecs.estimate, 'no'));
Settings.model.FlMode = [1 2 3];
Settings.model.objF = 'objf_LacI_TetR_simple';

data.dose = linspace(0, 2000000000, 1000)';
data.time = 1200;
data.tdh3 = 85.8657;
data.empty = 0.5198;

