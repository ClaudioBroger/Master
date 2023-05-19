clear all; clc;

%Reset, reload and compile model
sbioreset;
LacI_TetR_PMA1_model;
sbiosaveproject('LacI_TetRTup1_PMA1_model');
model = sbioloadproject('LacI_TetRTup1_PMA1_model');
num_draws = 30;

sbioaccelerate(model.mw_sbmod1)

%Modelsettings
ModelSettings_LacI_TetR_model_PMA1;

%draw parameter values
draw_parameters_LacI_TetR_PMA1;

%Simulate with drawn parameters
Simulate_random_drawn_parameters_TetRTup1_model;

%Drawing parameter values for EC50 ect calculations
num_draws = 200;
draw_parameters_LacI_TetR_with_aTc_EC50;

%Simulate with drawn parameters for EC50
Simulate_random_drawn_parametersTetRTup1_model_EC50;



rand_parameter.kTetR = randomValues;