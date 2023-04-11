%%  Sections to run
%%%%%%%%%%%%%%%%%%
draw_parameters             = 1;
plot_uncertainty            = 1;

model = sbioloadproject('LacImodel');
sbioaccelerate(model.mw_sbmod1)

ModelSettings_4;

p0_Fix = paramSpecs.p0(Settings.model.PIdxfixed);
for nFix=1:length(Settings.model.PIdxfixed)
    model.mw_sbmod1.Parameters(Settings.model.PIdxfixed(nFix)).Value = p0_Fix(nFix,1);
end

dataPos = strcat('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/code/results', "20_02_09_fitModel_Hyperspace_npoints20_Eline.mat");
load(dataPos)
num_draws = 10; %specify the number of random draws
%% draw_parameters

if draw_parameters
    paramSpecs = paramSpecs(Settings.model.PIdx,:);
    rand_parameter = draw_parameter_values(OutV,paramSpecs);
end

%% plot_uncertainty

if plot_uncertainty
    for k = 1:num_draws
        FlMode = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
        simulate_uncertainty_analysis(rand_parameter(:,k), model, Settings.model.dataPath, true, FlMode, Settings.model.PIdx, false, k);
    end
end



