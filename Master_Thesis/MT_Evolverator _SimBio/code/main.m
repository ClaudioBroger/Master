 clear all; clc;
% sbioreset
%% Sections to run
%%%%%%%%%%%%%%%%%
estimate_parameters         = 1;
select_optimal_param        = 1;
plot_meigo                  = 1;
runHyperspace				= 1;
plot_param_post             = 1;

%initPaths; to add dependencies (to do)

%% Load Model

% Step 1: create model in a seperate model file (e.g. TetR_model_dimensionless), run the model & save the model using  sbiosaveproject
% projectname (e.g.  sbiosaveproject TetRmodel)

% Step 2: Load model 
model = sbioloadproject('TetRmodel');

% Step 3: accelerate model
sbioaccelerate(model.mw_sbmod1)

% Step 4: Load Settings
ModelSettings; % note: you have to fill in the correct model, objF and parameterfiles in this settings file 

% Step 5: Set fixed parameters in model (according the parameter excel
% file).
p0_Fix = paramSpecs.p0(Settings.model.PIdxfixed);
for nFix=1:length(Settings.model.PIdxfixed)
    model.mw_sbmod1.Parameters(Settings.model.PIdxfixed(nFix)).Value = p0_Fix(nFix,1);
end

%% estimate parameters
if estimate_parameters
    %p0 = paramSpecs.p0(Settings.model.PIdx);
    p0 = ParaValues(Settings.model.PIdx);
    [~, ndata, ~, ~] = objf_TetR_experiments_0(p0, model, Settings.model.dataPath, false, Settings.model.FlMode, Settings.model.PIdx, false);
    threshold = 0.5*chi2inv(0.95,ndata - length(Settings.model.PIdx));
    [paraoptlist, costlist] = find_initial_points(Settings.model.objF, model, Settings.model.dataPath, false, Settings.model.FlMode, paramSpecs,Settings.model.PIdx, Settings.meigo.npoints, Settings.meigo.maxeval, threshold, false);
    savetofile(Settings.meigo.outputfile, paraoptlist, costlist, Settings);  
end
%% paraopt
% select optimal parameter set
if select_optimal_param
    %load(param_path);
    [mincost,idx_mincost] = min(costlist);
    paraopt = paraoptlist(idx_mincost,:);
end

%% plot optimal parameter set
if plot_meigo
    [objfc, ndata, objf, objfn] = objf_TetR_experiments_0(paraopt, model, Settings.model.dataPath, true, Settings.model.FlMode, Settings.model.PIdx, false);
end
%% uncertainty of parameters with hyperspace
if runHyperspace
    p0 = paramSpecs.p0(Settings.model.PIdx);
    [~, ndata, ~, ~] = objf_TetR_experiments_0(p0, model, Settings.model.dataPath, false, Settings.model.FlMode, Settings.model.PIdx, false);
    %threshold = chi2inv(0.95,ndata - length(Settings.model.PIdx));
    threshold = 3;
    OutV = sample_posterior(paraoptlist, Settings.model.objF, paramSpecs, Settings.model.PIdx, Settings.hyperspace.nsamples, threshold, model,  Settings.model.dataPath, false, Settings.model.FlMode, false);
    % save results: save(Settings.hyperspace.outputfile,'paraoptlist','costlist','Settings.model.FlMode','Settings.model.PIdx','OutV');
end

%% plot parameter posterior
if plot_param_post
    %load(hyperspace_load_path)
    figure;
    plot_OutV(OutV,paramSpecs(Settings.model.PIdx,:));
end

