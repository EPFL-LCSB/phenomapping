% Initialization of paths and loading and preparing the model for the
% phenomapping test
%
% USAGE:
%
%    run('settings.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

clear
clc

% Define a directory where we want to save final resutls (.mat and .csv) -
% if it does not exist it will be created (for this command to work you
% need to run this script / do not copy and paste!
% saving_directory = 'tmpresults/';
saving_directory = strrep(mfilename('fullpath'),...
    'tests/settings_itgo','tmpresults/');

% Initialize paths to phenomapping and dependencies (matTFA, texfba, and
% cplex)
addpath(genpath(strrep(mfilename('fullpath'),...
    'tests/settings_itgo','phenomapping')));
[phenomapping_directory, thermo_data_directory] = initPhenoMappingPaths(...
    saving_directory);

% Load test model: ipbe (tfa structure)
model = load(strcat(phenomapping_directory,'/models/itgo.mat'));
model = model.itgo;
modeldescription = 'iTgo';

% Prepare the model for phenomapping
[model, checkList] = initTestPhenoMappingModel(model,...
    thermo_data_directory);
clear thermo_data_directory

% Define path to phenotypes and omics data for the test model
phenotypes_description = 'phenotypes tachyzoites';
phenotypes_directory = which('itgo_tachy_phenotypes.mat');

metabolomics_description = {'metabolomics blood pfa'};
metabolomics_directory = which('allmetab_pfa_blood.mat');

transcriptomics_description = {'transcriptomics tachyzoites',...
    'transcriptomics bradyzoites'};
transcriptomics_directory{1} = which(...
    'itgo_tachy_levelGenes.mat');
transcriptomics_directory{2} = which(...
    'itgo_brady_levelGenes.mat');

filename = strcat(modeldescription,'_PhenoMappingSettings');
save(strcat(saving_directory,filename,'.mat'));
