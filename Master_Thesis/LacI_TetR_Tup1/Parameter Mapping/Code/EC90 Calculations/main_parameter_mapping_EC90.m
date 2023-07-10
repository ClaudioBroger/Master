clear all; clc;

%Reset, reload and compile model
sbioreset;
simple_model_LacI_TetR;
sbiosaveproject('simple_model');
model = sbioloadproject('simple_model');

sbioaccelerate(model.mw_sbmod1)

%Load parameters & SimfluoValues

load('27-Apr-2023_rand_parameter_LacI_TetR_EC50.mat');
load('27-Apr-2023SimFluoValues1_LacI_TetRTup1_EC50.mat');
load("27-Apr-2023SimFluoValues2_LacI_TetRTup1_EC50.mat");
load("27-Apr-2023SimFluoValues3_LacI_TetRTup1_EC50.mat");

%Modelsettings
ModelSettings_LacI_TetR_EC50;



%EC90 calculations
Ec90_calculations;

%Statistics EC90
statistics_ec90_TetRTup1;

%Order EC90
EC90_order_TetRTup1;

%All in the right order -> no statistics





