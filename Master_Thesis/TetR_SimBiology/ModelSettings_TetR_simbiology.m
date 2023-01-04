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

Settings.meigo.maxeval = 5000;
Settings.meigo.npoints = 1; %20 
Settings.meigo.outputfile =  ['Results' filesep '1102022_fitModel'];
%FlPlot = false;

% Hyperspace Settings

Settings.hyperspace.nsamples = 10000;
Settings.hyperspace.outputfile = 'Results\1102022_fitModel_hyperspace' ;

% Posterior
Settings.hyperspace.outputfile2 = 'Results\1102022_paramDist' ;