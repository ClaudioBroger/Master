%%% Make IQM model: Compiles the IQM model
function []= makeIQMmodel(modelfile, paramSpecs)

% Add IQMtools to path directory
%IQMdirectory=['..',filesep,'..',filesep,'..',filesep,'Dependencies',filesep,'IQMtools'];
% addpath(genpath(IQMdirectory))


[exdir,~,~]=fileparts(which(modelfile));

% change into right folder
OldFolder = cd;
cd(exdir);

% IQMmodel
model=IQMmodel(modelfile);

% Parameters
 % Set fixed parameters before compiling and estimated parameters
PIdxFixed = (strcmp(paramSpecs.estimate, 'no')); 
fixedParams = paramSpecs.names(strcmp(paramSpecs.estimate, 'no')); 
model = IQMparameters(model, fixedParams, paramSpecs.p0(PIdxFixed));

% Mex Model
IQMmakeMEXmodel(model, 'MEX_TetR');

% change back to old foler
cd(OldFolder);
end