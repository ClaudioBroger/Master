 clear all; clc;

%% Sections to run
%%%%%%%%%%%%%%%%%
create_model                = 1;
estimate_parameters         = 1;
select_optimal_param        = 1;
plot_meigo                  = 1;
runHyperspace				= 1;
plot_param_post             = 1;

ModelSettings;
%initPaths; to add dependencies

%% Create Model
if create_model
   disp('Making simbiology model')
   makeIQMmodel(Settings.model.modelfile, paramSpecs); % To do: Give mex name (from modelsettings) as input
end


%% estimate parameters
if estimate_parameters
    p0 = paramSpecs.p0(Settings.model.PIdx);
    [~, ndata, ~, ~] = objf_TetR_experiments_0(p0, Settings.model.fnHandle, Settings.model.dataPath, false, Settings.model.FlMode, Settings.model.PIdx, false);
    threshold = 0.5*chi2inv(0.95,ndata - length(Settings.model.PIdx));
    [paraoptlist, costlist] = find_initial_points(Settings.model.objF, Settings.model.fnHandle, Settings.model.dataPath, false, Settings.model.FlMode, paramSpecs,Settings.model.PIdx, Settings.meigo.npoints, Settings.meigo.maxeval, threshold, false);
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
    [objfc, ndata, objf, objfn] = objf_TetR_experiments_0(paraopt, Settings.model.fnHandle, Settings.model.dataPath, true, Settings.model.FlMode, Settings.model.PIdx, false);
end
%% uncertainty of parameters with hyperspace
if runHyperspace
    p0 = paramSpecs.p0(Settings.model.PIdx);
    [~, ndata, ~, ~] = objf_TetR_experiments_0(p0, Settings.model.fnHandle, Settings.model.dataPath, false, Settings.model.FlMode, Settings.model.PIdx, false);
    %threshold = chi2inv(0.95,ndata - length(Settings.model.PIdx));
    threshold = 1;
    OutV = sample_posterior(paraoptlist, Settings.model.objF, paramSpecs, Settings.model.PIdx, Settings.hyperspace.nsamples, threshold, Settings.model.fnHandle,  Settings.model.dataPath, false, Settings.model.FlMode, false);
    % save results: save(Settings.hyperspace.outputfile,'paraoptlist','costlist','Settings.model.FlMode','Settings.model.PIdx','OutV');
end

%% plot parameter posterior
if plot_param_post
    %load(hyperspace_load_path)
    figure;
    plot_OutV(OutV,paramSpecs(Settings.model.PIdx,:));
end

