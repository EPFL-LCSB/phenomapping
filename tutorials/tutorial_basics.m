%% Setting up the paths, solver and preparing the model for PhenoMapping
% files to set up the paths and prepare the model for PhenoMapping have
% been provided in the phenomapping/tests folder

% prior to running phenomapping you should install matTFA.
% TEX-FBA is also required should you want to integrate transcriptomics
% data and identify potential regulation of gene expression between
% isoenzymes

% it is recommended that you remove any other file from the matlab path 
% prior to starting phenomapping

% OPTION A: add the required repositories and solver directories to the 
% matlab path in an automatic fashion, meaning that phenomapping will 
% automatically search for matTFA and TEX-FBA in the parent folder and if 
% these are not there you will be asked to provide the paths. 
% For this option run the "settings.m" script

% OPTION B: add the required repositories and solver directories to the 
% matlab path (same as A) AND prepare the model for phenomapping analyses 
% - also done in an automatic fashion. 
% For this option run the "settings_nameofmodel.m" script

% OPTION C: add manually the required repositories and solver paths to the 
% matlab path AND prepare the model for phenomapping analyses. For this
% option you can follow the backbone of the code presented below. 
% Please, note that since option C is manual, you should adapt the paths

clear
clc
close all

% adding phenomapping and dependent repositories to the path
cplex_directory = '/Users/Anush/Applications/IBM/ILOG/CPLEX_Studio1271'; % provide path to cplex
phenomapping_directory = '/Users/Anush/GIT_Folders/phenomapping'; %provide path to phenomapping repository
mattfa_directory = '/Users/Anush/GIT_Folders/matTFA'; % provide path to matTFA repository
texfba_directory = '/Users/Anush/GIT_Folders/texfba'; % provide path to TEXFBA repository

addpath(genpath(phenomapping_directory));
addpath(genpath(mattfa_directory));
addpath(genpath(texfba_directory));
cd(phenomapping_directory) % work from first folder of phenomapping as the starting path

% providing information for TFA studies
tagThermo = 1;
thermo_data_directory = '/Users/Anush/GIT_Folders/matTFA/thermoDatabases/thermo_data.mat';

% providing the model
modeldescription = 'iPbe';
modelPath = '/Users/Anush/GIT_Folders/phenomapping/models/tipbe4liver.mat';
load(modelPath)
model = tipbe;

% path to save results
pathToSave = 'tmpresults/';

% preparing model for phenomapping
[model, checkList, tagReady] = initTestPhenoMappingModel(model,cplex_directory,thermo_data_directory);

save(strcat(pathToSave,'inputPhenoMapping.mat'), 'model');

%% PhenoMapping analysis for condition specific information in the GEM
% as applied for the analysis of iPbe and construction of iPbe-liver and 
% iPbe-blood 

% Step 0: Get essentiality for the non-stage specific or generic model. 
% This should be a model without context-specific data integrated. 
% iPbe is an example of such a generic model. 
% iPbe is saved in the files tipbe4blood.m and tipbe4liver.m. Both versions
% are stoichiometrically identical. They differ in the thermodynamic data
% integrated, i.e., they integrate pH and membrane potentials for the blood- 
% and liver-stages of the malaria infection. This allows that the 
% consequent phenomapping analyses are thermodynamically consistent with
% the corresponding life stage.
cd('tests')
run('test_core_essentiality.m')

% Step 1: Map essentiality to substrate availability following the 
% PhenoMapping analysis for the medium:
cd('tests')
run('test_core_substrates.m')

% Step 2: Map essentiality to metabolite concentrations following the 
% PhenoMapping analysis for metabolomics (if available):
cd('tests')
run('test_core_metabolomics.m')

% Step 3: Map essentiality to gene expression levels following the 
% PhenoMapping analysis for transcriptomics (if available) considering or 
% not regulation of gene expression between isoenzymes:
cd('tests')
run('test_core_mattfa_minmax.m') % this is to perform TVA on the model
% reaction fluxes as an optional input for TEX-FBA
run('test_core_transcriptomics.m')

%% Extract data stored in tmpresults
% The files stores in mat structures in tmpresults can be read an the 
% information extracted into text or csv files for further manual analysis
cd('tests')
run('test_io.m')