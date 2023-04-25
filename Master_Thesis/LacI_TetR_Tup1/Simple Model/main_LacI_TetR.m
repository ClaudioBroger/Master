clear all; clc;

%Reset, reload and compile model
sbioreset;
simple_model_LacI_TetR;
sbiosaveproject('simple_model');
model = sbioloadproject('simple_model');

sbioaccelerate(model.mw_sbmod1)

%Modelsettings
ModelSettings_LacI_TetR;

p0_Fix = paramSpecs.p0(Settings.model.PIdxfixed);
for nFix=1:length(Settings.model.PIdxfixed)
    model.mw_sbmod1.Parameters(Settings.model.PIdxfixed(nFix)).Value = p0_Fix(nFix,1);
end

p0 = ParaValues(Settings.model.PIdx);
[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
[paraValues] = cell2mat(get(model.mw_sbmod1.parameters, {'Value'}));
 paraValues = paramSpecs.p0;

 %Simulate model with p0
[objfc, ndata, objf, objfn, SimFluoValues1, SimFluoValues2, SimFluoValues3] = objf_LacI_TetR_simple(p0, model, Settings.model.PIdx, false, paraValues, ParaNames, data);
   
