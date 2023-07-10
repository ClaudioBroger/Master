clear all; clc;

%Reset, reload and compile model
sbioreset;
simple_model_LacI_TetR;
sbiosaveproject('simple_model');
model = sbioloadproject('simple_model');

sbioaccelerate(model.mw_sbmod1)

%Modelsettings
ModelSettings_LacI_TetR_EC10;

%draw parameter values
draw_parameters_LacI_TetR_EC50;

%Simulate with drawn parameters
Simulate_random_drawn_parameters_EC50;

%EC50 calculations
Ec10_calculations;

%Statistics EC10
statistics_ec10_TetRTup1;

%Order EC10
EC10_order_TetRTup1;

%All in the right order -> no statistics





