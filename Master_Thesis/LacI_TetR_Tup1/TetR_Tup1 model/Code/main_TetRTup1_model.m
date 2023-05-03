clear all; clc;

%Reset, reload and compile model
sbioreset;
LacI_TetRTup1_model_with_aTc;
sbiosaveproject('LacI_TetRTup1_model');
model = sbioloadproject('LacI_TetRTup1_model');

sbioaccelerate(model.mw_sbmod1)

%Modelsettings
ModelSettings_LacI_TetR_model;

%draw parameter values
draw_parameters_LacI_TetR_with_aTc;

%Simulate with drawn parameters
Simulate_random_drawn_parametersTetRTup1_model;

%Drawing parameter values for EC50 ect calculations
draw_parameters_LacI_TetR_with_aTc_EC50;

%Simulate with drawn parameters for EC50
Simulate_random_drawn_parametersTetRTup1_model_EC50;



