%%% Model & Estimation Settings %%%

Settings = [];

%General
Settings.model.modelproj = 'model.mw_sbmod1';
[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
[ParaValues] = cell2mat(get(model.mw_sbmod1.parameters, {'Value'}));
ParametersTable = readtable('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/parameters_0_LacI_TetR_model.xlsx');

% Specific
[~,Locb] = ismember(ParaNames, ParametersTable.names); %tells you how paramspecs should be changed
paramSpecs = ParametersTable(Locb, :);
Settings.model.PIdx = find(strcmp(paramSpecs.estimate, 'yes'));
Settings.model.PIdxfixed = find(strcmp(paramSpecs.estimate, 'no'));
Settings.model.FlMode = [1 2 3];
Settings.model.objF = 'objf_LacI_TetR_simple';

data_IPTG.dose = logspace(0,9,1000)';
data_IPTG.time = 1200;
data_IPTG.tdh3 = 85.8657;
data_IPTG.empty = 0.5198;

load("dataSRtup1.mat")
data.dose = data.dose*3.3*(1/480.9);
%data.dose = linspace(min(data.dose), max(data.dose), 1000);
