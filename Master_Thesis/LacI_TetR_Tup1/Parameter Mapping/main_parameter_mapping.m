clear all; clc;

%Reset, reload and compile model
sbioreset;
simple_model_LacI_TetR;
sbiosaveproject('simple_model');
model = sbioloadproject('simple_model');

sbioaccelerate(model.mw_sbmod1)

%Modelsettings
ModelSettings_LacI_TetR;

%draw parameter values
draw_parameters_LacI_TetR;

%Simulate with drawn parameters
Simulate_random_drawn_parameters;



