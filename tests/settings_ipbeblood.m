% Initialization of paths and loading and preparing the model for the
% phenomapping test
%
% USAGE:
%
%    run('settings_ipbeblood.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

clear
clc

% Define a directory where we want to save final resutls (.mat and .csv/
% .txt) - if the directory does not exist it will be created automatically
% Note: run this script (do not copy and paste in the command window)
saving_directory = strrep(mfilename('fullpath'),...
    'tests/settings_ipbeblood','tmpresults/');

% Initialize paths to phenomapping and dependencies (matTFA, texfba, and
% cplex)
addpath(genpath(strrep(mfilename('fullpath'),...
    'tests/settings_ipbeblood','phenomapping')));
[phenomapping_directory, thermo_data_directory] = initPhenoMappingPaths(...
    saving_directory);

% Load test model: ipbe (tfa structure)
model = load(strcat(phenomapping_directory,'/models/tipbe4blood.mat'));
model = model.tipbe;
modeldescription = 'iPbeBlood';

% Prepare the model for phenomapping
[model, checkList] = initTestPhenoMappingModel(model,...
    thermo_data_directory);
clear thermo_data_directory

% Define path to phenotypes and omics data for the test model
phenotypes_description = {'phenotypes blood pbe','phenotypes liver pbe'};
phenotypes_directory{1} = which('pbe_phenotypes_blood_Nov16.mat');
phenotypes_directory{2} = which('pbe_phenotypes_liver_Jul18.mat');

metabolomics_description = {'metabolomics blood pfa'};
metabolomics_directory = which('allmetab_pfa_blood.mat');

transcriptomics_description = {'transcriptomics blood pbe'};
transcriptomics_directory = which(...
    'levelGenes_pbe_blood_Troph_mean_24h.mat');

filename = strcat(modeldescription,'_PhenoMappingSettings');
save(strcat(saving_directory,filename,'.mat'));
