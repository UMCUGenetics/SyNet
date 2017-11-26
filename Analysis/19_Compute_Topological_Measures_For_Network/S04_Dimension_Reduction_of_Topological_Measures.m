clc;
clear;

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
n_pair = 10000;
IS_STRICT_CV = 1;
CLS_Name = 'Lasso';

%% Load TM data
load('./Topological_Data/TMData_NS20000_NF90.mat', 'zTM_Data', 'TM_Label', 'TM_Name');

%%