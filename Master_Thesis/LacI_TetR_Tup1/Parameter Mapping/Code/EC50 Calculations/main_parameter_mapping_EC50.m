clear all; clc;

%Reset, reload and compile model
sbioreset;
simple_model_LacI_TetR;
sbiosaveproject('simple_model');
model = sbioloadproject('simple_model');

sbioaccelerate(model.mw_sbmod1)

%Modelsettings
ModelSettings_LacI_TetR_EC50;

%draw parameter values
draw_parameters_LacI_TetR_EC50;

%Simulate with drawn parameters
Simulate_random_drawn_parameters_EC50;

%EC50 calculations
Ec50_calculations;

%Statistics EC50
statistics_ec50_TetRTup1;

%Order EC50
EC50_order_TetRTup1;

%All in the right order -> no statistics





