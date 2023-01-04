%%%% Model & Estimation Settings %%%%%%%%%%%%%%%

Settings = [];

% General
Settings.model.modelfile = 'Modelv001.txt';
model = IQMmodel(Settings.model.modelfile); 
[ParaNames, ParaValues] = IQMparameters(model);
ParametersTable = readtable('parameters_0.xlsx');
Settings.model.dataPath = ('DataForEstimation\');

% Specific
[~,Locb] = ismember(ParaNames, ParametersTable.names); %tells you how paramspecs should be changed
paramSpecs = ParametersTable(Locb, :);
Settings.model.PIdx = find(strcmp(paramSpecs.estimate, 'yes')) ;
Settings.model.FlMode = [1];
Settings.model.fnHandle = @MEX_TetR;
Settings.model.objF = 'objf_TetR_experiments_0' ;

% Meigo Estimation Settings

Settings.meigo.maxeval = 3000;
Settings.meigo.npoints = 1; %20 
Settings.meigo.outputfile =  ['\\mac\home\Documents\ETH\Master_Thesis\Evolverator\MT_Evolverator\code\Results' filesep '1102022_fitModel'];
%FlPlot = false;

% Hyperspace Settings

Settings.hyperspace.nsamples = 3000;
Settings.hyperspace.outputfile = '\\mac\home\Documents\ETH\Master_Thesis\Evolverator\MT_Evolverator\code\Results\1102022_fitModel_hyperspace' ;

% Posterior
Settings.hyperspace.outputfile2 = '\\mac\home\Documents\ETH\Master_Thesis\Evolverator\MT_Evolverator\code\Results\1102022_paramDist' ;