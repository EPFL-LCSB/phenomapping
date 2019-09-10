%% Setting up the paths, solver and preparing the model for PhenoMapping
% files to set up the paths and prepare the model for PhenoMapping have
% been provided in the phenomapping/tests folder

% prior to running phenomapping you should install matTFA.
% TEX-FBA is also required should you want to integrate transcriptomics
% data and identify potential regulation of gene expression between
% isoenzymes

% Work from the parent phenomapping folder.

% it is recommended that you remove any other file from the matlab path 
% prior to starting phenomapping

% OPTION A: add the required repositories and solver directories to the 
% matlab path in an automatic fashion, meaning that phenomapping will 
% automatically search for matTFA and TEX-FBA in the parent folder and if 
% these are not there you will be asked to provide the paths. 
% For this option run the "settings.m" script.
% You should work from the parent phenomapping folder to run "settings.m"

% OPTION B: add the required repositories and solver directories to the 
% matlab path (same as A) AND prepare the model for phenomapping analyses 
% - also done in an automatic fashion. 
% For this option run the "settings_nameofmodel.m" script
% You should work from the parent phenomapping folder to run 
% "settings_nameofmodel.m"

% OPTION C: add manually the required repositories and solver paths to the 
% matlab path AND prepare the model for phenomapping analyses. For this
% option you can follow the backbone of the code presented below. 
% Please, note that since option C is manual, you should adapt the paths
% You should work from the parent phenomapping folder to run 
% "initTestPhenoMappingModel.m"

clear
clc
close all


% Provide path to folder where your results will be saved - please, keep 
% the name "tmpresults" for this folder.
% initPhenoMappingPaths will create the "tmpresults" folder at the same 
% level as ext, models, phenomapping, tests, and tutorials subfolders.
% Ideally you will provide here the full path to the saving directory.
% Something like:
saving_directory = '/Users/Anush/GIT_Folders/phenomapping/tmpresults/';
% Alternatively, if you work from the parent phenomapping folder you can
% define:
% saving_directory = 'tmpresults/';


% See in the comment below an example of the paths to phenomapping, matTFA
% and TEX-FBA repositories and CPLEX
phenomapping_directory = '/Users/Anush/GIT_Folders/phenomapping';
% mattfa_directory = '/Users/Anush/GIT_Folders/matTFA';
% texfba_directory = '/Users/Anush/GIT_Folders/texfba';
% cplex_directory = '/Users/Anush/Applications/IBM/ILOG/CPLEX_Studio1271';

addpath(genpath(phenomapping_directory));

% Call the function "initPhenoMappingPaths.m" to add all other required 
% paths
[phenomapping_directory, thermo_data_directory] = initPhenoMappingPaths(...
    saving_directory);

% NOTE that matTFA and TEX-FBA are in the same directory as phenomapping.
% This will be assumed by default in the "initPhenoMappingPaths.m".
% If this is not the case you will have to provide the path to matTFA and 
% TEX-FBA when you are asked to do so by "initPhenoMappingPaths.m".

% You will also need to paste the path to cplex when the following message
% appears:
% "Please provide your cplex path and press enter"
% You will paste the corresponding path as follows (without " "):
% ... /Users/Anush/Applications/IBM/ILOG/CPLEX_Studio1271


% Work from parent folder of phenomapping as the starting path. If you are
% not there run the following command:
% cd(phenomapping_directory)


% Tag true to apply phenomapping with TFA (thermodynamic constraints will
% be taken into account)
tagThermo = 1;

% Upload model - here iPbe for liver-stage analysis provided in the
% subfolder models as test case
modeldescription = 'iPbe';
modelPath = '/Users/Anush/GIT_Folders/phenomapping/models/tipbe4liver.mat';
load(modelPath)
model = tipbe;

% Prepare model for phenomapping
[model, checkList, tagReady] = initTestPhenoMappingModel(model,...
    thermo_data_directory,tagThermo);

% Save prepared model
save(strcat(saving_directory,'inputPhenoMapping.mat'), 'model');

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
