%%%% Model & Estimation Settings %%%%%%%%%%%%%%%

Settings = [];

% General
Settings.model.modelproj = 'model.mw_sbmod1';
[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
[ParaValues] = cell2mat(get(model.mw_sbmod1.parameters, {'Value'}));
ParametersTable = readtable('/model/parameters_0.xlsx');
Settings.model.dataPath = ('/Users/claudiobroger/Documents/ETH/Master_Thesis/MT_Evolverator _SimBio/code/DataForEstimation/');

% Specific
[~,Locb] = ismember(ParaNames, ParametersTable.names); %tells you how paramspecs should be changed
paramSpecs = ParametersTable(Locb, :);
Settings.model.PIdx = find(strcmp(paramSpecs.estimate, 'yes')) ;
Settings.model.PIdxfixed = find(strcmp(paramSpecs.estimate, 'no')) ;
Settings.model.FlMode = [1];
Settings.model.objF = 'objf_TetR_experiments_0' ;

% Meigo Estimation Settings

Settings.meigo.maxeval = 5000;
Settings.meigo.npoints = 1; %20 
Settings.meigo.outputfile =  ['/Users/claudiobroger/Documents/ETH/Master_Thesis/MT_Evolverator _SimBio/code/results' filesep '20122022_fitModel'];
%FlPlot = false;

% Hyperspace Settings

Settings.hyperspace.nsamples = 2000;
Settings.hyperspace.outputfile = ['/Users/claudiobroger/Documents/ETH/Master_Thesis/MT_Evolverator _SimBio/code/results' filesep '20122022_fitModel_hyperspace'] ;

% Posterior
Settings.hyperspace.outputfile2 = ['/Users/claudiobroger/Documents/ETH/Master_Thesis/MT_Evolverator _SimBio/code/results' filesep '20122022_paramDist'] ;