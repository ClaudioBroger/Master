%%% ----- Main script from which you can run everything --------- %%%%

%% Create Model
    % If desired, give model a specific name

model = IQMmodel('code/model/model.txtbc');
if ~isfile('code/model/model.mexw64')
IQMmakeMEXmodel(model, 'model');
end

%% Extract current parametervalues / names

[names, values] = IQMparameters(model) ;
ParaNames = names ;
para = values ; 

%% Run IQM model

tspan = [0:1:1400] ;
output = IQMPsimulate('model', tspan , [], ParaNames, para) ; 
        % [] represents initial conditions. Use [] when you want to keep IC
        % specified in the model file

% Quick plot to see simulation result
figure(1)
plot(output.time, output.statevalues)
figure(2)
plot(output.time, output.variablevalues)

%% Run objective function

para_est =  ;
FlMode = [1] ; % which module
FlPlot = true; % true = plot, false = no plot
PIdx_full = ;
%para(PIdx_full) = para_est ;

[objfc, ndata, objf, objfn] = objf_function(para_est, FlPlot, FlMode, PIdx_full);

%% Run optimization